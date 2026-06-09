% ========================================================================
%  FUNCIÓN: SolenoidFields.m
%
%  Propósito:
%    Calcula el campo magnético 3D generado por un solenoide de NL vueltas
%    usando la Ley de Biot-Savart discretizada sobre segmentos de corriente.
%    También calcula el campo eléctrico inducido en el plano transversal XY.
%
%  Entradas:
%    nl  - Número de vueltas del solenoide
%    ds  - Paso espacial de la grilla (m); menor = mayor resolución
%    rw  - Radio del cable conductor (m); evita singularidad en r=0
%    I   - Corriente eléctrica (Amperios)
%    N   - Puntos de discretización por espira
%    R   - Radio del solenoide (m)
%    sz  - Separación axial entre vueltas consecutivas (m)
%
%  Salidas:
%    Bx, By, Bz  - Componentes del campo magnético (matrices 3D) [Teslas]
%    Ex, Ey      - Componentes del campo eléctrico inducido (matriz 2D) [V/m]
%    Emag        - Magnitud del campo eléctrico (matriz 2D) [V/m]
%    x, y, z     - Vectores de coordenadas de la grilla [m]
%    Px, Py, Pz  - Coordenadas geométricas de los nodos del cable [m]
% ========================================================================
function [Bx, By, Bz, Ex, Ey, Emag, x, y, z, Px, Py, Pz] = SolenoidFields(nl, ds, rw, I, N, R, sz)

% -----------------------------------------------------------------------
%  BLOQUE 0: PREPARACIÓN DEL ESPACIO DE SIMULACIÓN
% -----------------------------------------------------------------------
% Se define una grilla cúbica de -5 a 5 metros con paso ds.
% ndgrid (usado más adelante) requiere vectores 1D consistentes en los 3 ejes.
x = (-5:ds:5); 
y = (-5:ds:5); 
z = x;                       % Grilla cúbica: mismo rango en Z que en X

Lx = length(x); Ly = length(y); Lz = length(z);

mo = 4*pi*1e-7;              % Permeabilidad magnética del vacío (H/m)
km = mo*I/(4*pi);            % Constante de Biot-Savart: (μ₀·I)/(4π)
                             % Agrupa la corriente y la constante para
                             % simplificar el cálculo en el bucle principal

% Encontrar el índice del plano Z=0 (plano central transversal).
% Se usa para extraer el campo eléctrico inducido en ese corte.
cp = find(z<=0); 
cz = cp(1);                  % Primer índice donde z <= 0 (más cercano a z=0)

% Parámetros de discretización angular del embobinado
s      = 1;                  % Puntero de índice para llenar Px, Py, Pz
dtheta = 2*pi/N;             % Ángulo entre nodos consecutivos de una espira
ang    = (0:dtheta:2*pi-dtheta); % Vector de ángulos para una espira completa (N puntos)

% -----------------------------------------------------------------------
%  BLOQUE 1: CONSTRUCCIÓN GEOMÉTRICA DEL SOLENOIDE
%
%  Se discretiza el cable del solenoide en nl*N segmentos rectilíneos.
%  Cada segmento queda definido por:
%    - Su punto de origen (Px, Py, Pz)
%    - Su vector diferencial de longitud (dx, dy, dz) = dL
%
%  Este vector dL es el elemento fundamental de la Ley de Biot-Savart:
%    dB = (μ₀I/4π) * (dL × r̂) / |r|²
% -----------------------------------------------------------------------
Px = zeros(1, N*nl); Py = zeros(1, N*nl); Pz = zeros(1, N*nl);
dx = zeros(1, N*nl); dy = zeros(1, N*nl); dz = zeros(1, N*nl);

for i_loop = 1:nl
    idx = s:s+N-1;           % Índices correspondientes a la vuelta i_loop

    % Posición de cada nodo sobre la circunferencia de radio R
    Px(idx) = R*cos(ang);
    Py(idx) = R*sin(ang);

    % Posición axial base de la vuelta (centrada en z=0)
    Pz(idx) = -nl/2*sz + (i_loop-1)*sz;

    % Diferencial de longitud en el plano XY (tangente a la circunferencia)
    % Deriva de la parametrización: d/dθ[R·cos(θ)] = -R·sin(θ) = -Py
    %                               d/dθ[R·sin(θ)] =  R·cos(θ) =  Px
    dx(idx) = -Py(idx)*dtheta;
    dy(idx) =  Px(idx)*dtheta;

    % Para nl > 1: distribuir los nodos a lo largo del paso axial sz
    % Esto convierte espiras planas en una hélice continua real.
    % Sin esto, todos los nodos de una vuelta tendrían la misma Z,
    % ignorando el avance axial del embobinado.
    if nl > 1
        Pz(idx) = Pz(idx) + linspace(0, N-1, N)*sz/N;
    end

    s = s + N;               % Avanzar el puntero a la siguiente vuelta
end

% Componente axial del diferencial de longitud dz
% Para una sola espira (nl=1): el cable es plano, no hay avance en Z
% Para múltiples espiras: cada segmento avanza sz/N metros en Z
if nl == 1
    dz(:) = 0;
else
    dz(:) = sz/N;
end

% -----------------------------------------------------------------------
%  BLOQUE 2: CÁLCULO DEL CAMPO MAGNÉTICO 3D — LEY DE BIOT-SAVART
%
%  La Ley de Biot-Savart para un segmento de corriente es:
%
%    dB = (μ₀I/4π) * (dL × r) / |r|³
%
%  donde r = punto_campo - punto_fuente es el vector de separación.
%
%  Implementación vectorizada:
%    - ndgrid genera matrices 3D [Lx × Ly × Lz] para X, Y, Z
%    - El bucle itera solo sobre los nl*N segmentos del cable
%    - Para cada segmento, rx/ry/rz son matrices 3D completas,
%      permitiendo calcular la contribución en todo el espacio de una vez
%    - Complejidad: O(nl·N·Lx·Ly·Lz) en lugar de O(nl·N·Lx·Ly·Lz) con bucles anidados
%      (misma complejidad, pero MATLAB ejecuta las operaciones matriciales en C)
%
%  Regularización: se suma rw² al denominador para evitar la singularidad
%  física en r=0 (sobre el cable). rw representa el radio finito del conductor.
% -----------------------------------------------------------------------
[X, Y, Z] = ndgrid(x, y, z);   % Grilla 3D con convención [X, Y, Z] estricta
Bx = zeros(Lx, Ly, Lz); 
By = zeros(Lx, Ly, Lz); 
Bz = zeros(Lx, Ly, Lz);

for L = 1:(nl*N)
    % Vector r = punto_campo - punto_fuente (para cada punto de la grilla)
    rx = X - Px(L);
    ry = Y - Py(L);
    rz = Z - Pz(L);

    % |r|³ regularizado: rw² evita la divergencia cuando r → 0
    r3 = (rx.^2 + ry.^2 + rz.^2 + rw^2).^(1.5);

    % Producto cruz dL × r, componente a componente:
    %   (dL × r)_x = dy·rz - dz·ry
    %   (dL × r)_y = dz·rx - dx·rz
    %   (dL × r)_z = dx·ry - dy·rx
    Bx = Bx + km * (dy(L).*rz - dz(L).*ry) ./ r3;
    By = By + km * (dz(L).*rx - dx(L).*rz) ./ r3;
    Bz = Bz + km * (dx(L).*ry - dy(L).*rx) ./ r3;
end

% -----------------------------------------------------------------------
%  BLOQUE 3: CAMPO ELÉCTRICO INDUCIDO EN EL PLANO TRANSVERSAL XY
%
%  En el plano central (z=0), el campo magnético axial Bz induce un campo
%  eléctrico azimutal por la Ley de Faraday. Para una distribución
%  cilíndrica, este campo eléctrico circula en sentido tangencial:
%
%    E · 2πr = -dΦ/dt   →   E_φ = -Bz · (dΦ/dt) / (2πr)
%
%  En coordenadas cartesianas, las componentes azimutales son:
%    Ex =  Bz · y / (2π·r)  ·  (factor de variación temporal)
%    Ey = -Bz · x / (2π·r)  ·  (factor de variación temporal)
%
%  Nota: el factor dtheta/2π es una aproximación proporcional al cambio
%  angular discretizado, no una derivada temporal explícita.
%  Esta formulación es válida cualitativamente para visualización.
% -----------------------------------------------------------------------

% Extraer el campo Bz en el plano Z central (corte transversal)
Bz_xy = squeeze(Bz(:,:,cz));   % Resultado: matriz 2D [Lx × Ly]

% Grilla 2D para el plano XY
[X2D, Y2D] = ndgrid(x, y);

% Radio en el plano transversal, regularizado con rw para evitar r=0
r_elec = sqrt(X2D.^2 + Y2D.^2 + rw^2);

% Componentes del campo eléctrico inducido (azimutal → cartesiano)
Ex   =  Bz_xy .* Y2D .* dtheta ./ (2*pi*r_elec);
Ey   = -Bz_xy .* X2D .* dtheta ./ (2*pi*r_elec);
Emag = sqrt(Ex.^2 + Ey.^2);    % Magnitud del campo eléctrico

end