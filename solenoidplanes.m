% ========================================================================
%  ARCHIVO CENTRALIZADOR: solenoidplanes.m
%
%  Propósito:
%    Genera todas las visualizaciones del campo magnético y eléctrico
%    producido por un solenoide, comparando 4 configuraciones de vueltas:
%    NL = 1, 3, 5 y 15.
%
%    Para cada configuración ejecuta 4 estudios:
%      Estudio 1 — Perfil longitudinal del campo B (plano XZ)
%      Estudio 2 — Campo eléctrico inducido (plano XY)
%      Estudio 3 — Perfiles 1D de componentes de B y E
%      Estudio 4 — Superficie 3D de la magnitud |B|
%
%    Al finalizar el bucle, genera el Estudio 5:
%      Estudio 5 — Geometría 3D real del embobinado con flechas de corriente
%
%  Dependencia:
%    Requiere la función SolenoidFields.m en el mismo directorio.
% ========================================================================
clear; clc; clf;

fprintf('=== INICIANDO VISUALIZACIÓN CENTRALIZADA DE CAMPOS Y GEOMETRÍA ===\n\n');

% -----------------------------------------------------------------------
%  PARÁMETROS GLOBALES DEL SOLENOIDE
%  Todos los estudios usan estos mismos valores para garantizar
%  consistencia física entre las comparaciones.
% -----------------------------------------------------------------------
NL_vector = [1, 3, 5, 15]; % Casos de estudio: número de vueltas a comparar
ds = 0.15;                  % Paso espacial de la grilla (m); 0.15 = buena resolución
rw = 0.2;                   % Radio del cable (m); regulariza la singularidad en r=0
I  = 3;                     % Corriente eléctrica en el conductor (Amperios)
N  = 20;                    % Puntos de discretización por espira (para Biot-Savart)
R  = 1.5;                   % Radio del solenoide (m)
sz = 0.35;                  % Separación axial entre vueltas consecutivas (m)


% ========================================================================
%  BUCLE PRINCIPAL: itera sobre cada configuración de NL
% ========================================================================
for H = 1:4
    nl_actual = NL_vector(H);
    fprintf('Procesando campo magnético para NL = %d...\n', nl_actual);

    % Llamada a SolenoidFields: calcula el campo B y E en toda la grilla 3D
    % usando la Ley de Biot-Savart sobre nl_actual*N segmentos de corriente.
    % dBx, dBy, dBz : componentes del campo magnético [Teslas] (matrices 3D)
    % Ex, Ey, Emag  : campo eléctrico inducido en plano XY (matrices 2D)
    % x, y, z       : vectores de coordenadas de la grilla
    % Px, Py, Pz    : geometría del cable del solenoide
    [dBx, dBy, dBz, Ex, Ey, Emag, x, y, z, Px, Py, Pz] = ...
        SolenoidFields(nl_actual, ds, rw, I, N, R, sz);

    % -------------------------------------------------------------------
    %  ESTUDIO 1: PERFIL LONGITUDINAL DEL CAMPO B — PLANO XZ
    %
    %  Se extrae el corte central en Y para obtener una vista de perfil
    %  del solenoide. Este plano muestra:
    %    - La distribución de |B| como mapa de calor (escala cúbica)
    %    - Las líneas de flujo magnético (streamslice)
    %    - La silueta del cable proyectada sobre el plano XZ
    %
    %  La escala cúbica (Bxz^(1/3)) comprime el rango dinámico para que
    %  las zonas de campo débil sean visibles junto a las de campo fuerte.
    % -------------------------------------------------------------------
    figure
    hold on

    % Magnitud total del campo en toda la grilla 3D
    Bmag = sqrt(dBx.^2 + dBy.^2 + dBz.^2);

    % Extraer el corte central en Y (plano de simetría del solenoide)
    centery = round(length(y)/2);
    Bx_xz = squeeze(dBx(:, centery, :)); % Componente radial en plano XZ
    Bz_xz = squeeze(dBz(:, centery, :)); % Componente axial en plano XZ
    Bxz   = squeeze(Bmag(:, centery, :)); % Magnitud en plano XZ

    % Mapa de calor: transponer para alinear Z vertical y X horizontal
    % La raíz cúbica comprime la escala sin perder información cualitativa
    pcolor(x, z, (Bxz').^(1/3));
    shading interp;
    colormap(jet);

    % Líneas de flujo del campo magnético en el plano XZ
    % Densidad 2 = balance entre legibilidad y detalle
    h1 = streamslice(x, z, Bx_xz', Bz_xz', 2);
    set(h1, 'Color', [0.9 1 0.9], 'LineWidth', 1.2);

    % Silueta del cable proyectada en XZ (línea discontinua negra)
    % Muestra dónde físicamente está el conductor respecto al campo
    dtheta    = 2*pi/N;
    dx_vector = -Py * dtheta;
    dz_vector = zeros(size(Pz));
    plot(Px, Pz, '--k', 'LineWidth', 1);

    axis([-5 5 -5 5]);
    grid on;
    title(['Campo B_{xz} (Perfil) para NL = ' num2str(nl_actual)]);
    xlabel('Eje X (Radio)'); ylabel('Eje Z (Longitud)');
    colorbar('Position', [0.93 0.11 0.02 0.81]);

    % -------------------------------------------------------------------
    %  ESTUDIO 2: CAMPO ELÉCTRICO INDUCIDO — PLANO XY
    %
    %  El campo magnético variable en el tiempo induce un campo eléctrico
    %  azimutal por la Ley de Faraday. Este estudio lo visualiza en el
    %  plano transversal central (Z=0).
    %
    %  El campo eléctrico calculado en SolenoidFields es una aproximación
    %  proporcional a Bz; útil para comparación cualitativa entre casos.
    % -------------------------------------------------------------------
    fprintf('\nGenerando mapa de Campo Eléctrico para NL = %d...\n', nl_actual);
    figure
    pcolor(x, y, Emag);
    shading interp;
    colormap(hot);
    colorbar;
    hold on

    % Líneas de flujo circulares del campo eléctrico inducido
    % El campo eléctrico es azimutal: sus líneas son circunferencias
    h2 = streamslice(x, y, Ex', Ey', 1.5);
    set(h2, 'Color', [0 1 0], 'LineWidth', 1.2);

    % Silueta frontal del solenoide (circunferencia de radio R)
    ang = linspace(0, 2*pi, 100);
    plot(R*cos(ang), R*sin(ang), 'w--', 'LineWidth', 2);

    axis([-R-2 R+2 -R-2 R+2]);
    axis square; grid on;
    title(['Campo Eléctrico Inducido E_{xy} (NL = ' num2str(nl_actual) ')']);
    xlabel('Eje X'); ylabel('Eje Y');

    % -------------------------------------------------------------------
    %  ESTUDIO 3: PERFILES 1D DE COMPONENTES DE B Y E
    %
    %  Las gráficas de línea permiten análisis cuantitativo que los mapas
    %  de calor no pueden dar:
    %
    %  Subplot 1 — Bz(z) en el eje central: mide la uniformidad del campo
    %    axial. Un plateau plano indica buen confinamiento longitudinal.
    %    Las líneas rojas marcan los extremos físicos del solenoide.
    %
    %  Subplot 2 — Bx(x) en el borde axial: muestra el campo radial en el
    %    extremo del solenoide, donde el campo pierde confinamiento.
    %    Un Bx grande indica fuga axial del campo.
    %
    %  Subplot 3 — Ey(x): componente del campo eléctrico inducido en Y,
    %    evaluada a lo largo del eje X. Máximo esperado en r = R.
    %
    %  Subplot 4 — Ex(y): componente del campo eléctrico en X,
    %    evaluada a lo largo del eje Y. Antisimétrica respecto a y=0.
    % -------------------------------------------------------------------
    fprintf('\nGenerando perfiles 1D de las componentes físicas...\n');
    figure

    % Índices del punto central de la grilla (x=0, y=0, z=0)
    cx = round(length(x)/2);
    cy = round(length(y)/2);
    cz = round(length(z)/2);

    % --- Subplot 1: Bz a lo largo del eje Z ---
    subplot(2,2,1)
    Bz_z = squeeze(dBz(cx, cy, :));
    plot(z, Bz_z, 'b-', 'LineWidth', 2);
    grid on; hold on;
    L_mitad = (nl_actual * sz) / 2;
    xline(-L_mitad, 'r--', 'Inicio');
    xline( L_mitad, 'r--', 'Fin');
    title('Campo Longitudinal (B_z) a lo largo de Z');
    xlabel('Eje Z (m)'); ylabel('B_z (T)');

    % --- Subplot 2: Bx a lo largo de X en el borde axial ---
    subplot(2,2,2)
    z_borde_idx = find(z >= L_mitad, 1); % Índice del extremo del solenoide
    Bx_x = squeeze(dBx(:, cy, z_borde_idx));
    plot(x, Bx_x, 'r-', 'LineWidth', 2);
    grid on; hold on;
    xline(-R, 'k:', 'Pared'); xline(R, 'k:', 'Pared');
    title(['Campo Radial (B_x) vs X (Corte en Z = ' num2str(z(z_borde_idx)) 'm)']);
    xlabel('Eje X (m)'); ylabel('B_x (T)');

    % --- Subplot 3: Ey a lo largo de X ---
    subplot(2,2,3)
    Ey_x = squeeze(Ey(:, cy));
    plot(x, Ey_x, 'g-', 'LineWidth', 2);
    grid on; hold on;
    xline(-R, 'k:', 'Pared'); xline(R, 'k:', 'Pared');
    title('Componente Eléctrica E_y a lo largo del Eje X');
    xlabel('Eje X (m)'); ylabel('E_y (V/m)');

    % --- Subplot 4: Ex a lo largo de Y ---
    subplot(2,2,4)
    Ex_y = squeeze(Ex(cx, :));
    plot(y, Ex_y, 'm-', 'LineWidth', 2);
    grid on; hold on;
    xline(-R, 'k:', 'Pared'); xline(R, 'k:', 'Pared');
    title('Componente Eléctrica E_x a lo largo del Eje Y');
    xlabel('Eje Y (m)'); ylabel('E_x (V/m)');

    % -------------------------------------------------------------------
    %  ESTUDIO 4: SUPERFICIE 3D DE LA MAGNITUD |B|
    %
    %  Visualiza la topología completa del campo magnético como un paisaje
    %  tridimensional donde la altura representa la intensidad de |B|.
    %
    %  Interpretación física:
    %    - La "meseta" central alta corresponde al interior del solenoide,
    %      donde el campo es uniforme y fuerte → zona de confinamiento.
    %    - Las "laderas" que caen hacia los bordes muestran dónde el campo
    %      se debilita → zona de fuga potencial.
    %    - A mayor NL, la meseta es más plana y ancha: mejor confinamiento.
    % -------------------------------------------------------------------
    figure
    cy = round(length(y)/2);

    % Magnitud del campo en el plano XZ central
    Bmag_xz = squeeze(sqrt(dBx(:,cy,:).^2 + dBy(:,cy,:).^2 + dBz(:,cy,:).^2));

    % surf con transparencia parcial para apreciar la estructura interior
    surf(z, x, Bmag_xz, 'FaceAlpha', 0.75);
    shading interp;
    colormap(jet);
    colorbar;

    title('Topología 3D de la Magnitud del Campo Magnético |B|');
    xlabel('Eje Z (Longitud del Solenoide en m)');
    ylabel('Eje X (Radio del Solenoide en m)');
    zlabel('Magnitud |B| (Teslas)');
    view(45, 35);
    axis tight;
    grid on;

end % Fin del bucle sobre NL_vector


% ========================================================================
%  ESTUDIO 5: GEOMETRÍA 3D DEL EMBOBINADO — TRAYECTORIA REAL DEL CABLE
%
%  Visualiza la forma física del solenoide para cada NL en subplots
%  comparativos. Muestra tres elementos superpuestos:
%    1. Línea negra continua: trayectoria helicoidal del cable
%    2. Puntos azules: nodos donde se aplica Biot-Savart
%    3. Flechas rojas: dirección y sentido de la corriente en cada segmento
%
%  Esta figura es útil como figura introductoria en el reporte para
%  contextualizar físicamente de dónde proviene el campo calculado.
%  La escala es idéntica en los 4 subplots para comparación directa.
% ========================================================================
fprintf('\nGenerando visualización tridimensional de la trayectoria de los loops...\n');
figure

for H = 1:4
    nl_actual = NL_vector(H);

    % Solo se necesita la geometría del cable; los campos se descartan (~)
    [~, ~, ~, ~, ~, ~, ~, ~, ~, Px, Py, Pz] = ...
        SolenoidFields(nl_actual, ds, rw, I, N, R, sz);

    % Vectores diferenciales de dirección de corriente en cada segmento
    % dx, dy: componentes tangenciales en el plano XY (circulación)
    % dz    : componente axial (avance helicoidal entre vueltas)
    dtheta    = 2*pi/N;
    dx_vector = -Py * dtheta;
    dy_vector =  Px * dtheta;
    dz_vector = zeros(size(Pz));
    if nl_actual > 1
        dz_vector(:) = sz/N; % Avance axial por segmento
    end

    subplot(1, 4, H)
    hold on

    % Trayectoria continua del cable
    plot3(Px, Py, Pz, '-k', 'LineWidth', 1.5);

    % Nodos de discretización (donde se evaluó Biot-Savart)
    plot3(Px, Py, Pz, '.b', 'MarkerSize', 8);

    % Flechas de dirección de corriente en cada segmento
    quiver3(Px, Py, Pz, dx_vector, dy_vector, dz_vector, 0.8, 'r', 'LineWidth', 1.2);

    grid on;
    view(-35, 40);
    axis([-R-0.5 R+0.5 -R-0.5 R+0.5 -4 4]);
    daspect([1 1 1]); % Proporción 1:1:1 para espiras perfectamente circulares

    title(['Loops para NL = ' num2str(nl_actual)]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
end

fprintf('\n=== TODAS LAS GRÁFICAS GENERADAS CON ÉXITO ===\n');