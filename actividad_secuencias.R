###############################################################################
#
#   ACTIVIDAD: DEL GEN A LA PROTEÍNA EN VARIAS ESPECIES
#   Genética Molecular Evolutiva Humana · Universidad del Tolima
#
###############################################################################
#
#   ¿QUÉ VAS A HACER?
#   1. Descargar del NCBI el MISMO gen en 3 a 5 especies (genes ortólogos).
#   2. Alinear las secuencias y ver con colores qué regiones se conservan.
#   3. Identificar el promotor, el inicio y el fin de la transcripción,
#      los exones, los intrones y los límites exón-intrón.
#   4. Transcribir un exón (ADN -> ARNm).
#   5. Traducir ese exón (ARNm -> proteína).
#   6. Comparar ADN, ARNm y proteína en un tramo donde las especies difieren.
#   7. Introducir una mutación y ver su efecto en la proteína.
#   8. Ver la proteína en 3D, con sus aminoácidos y mutaciones marcados.
#   9. Responder preguntas, hacer búsquedas y revisar tus respuestas.
#
#   INSTRUCCIONES
#   • Guarda ESTE archivo y "funciones_secuencias.R" en la MISMA carpeta.
#   • Solo modifica las secciones marcadas con  ✏️ .
#   • Ejecuta el script por secciones: selecciona las líneas y presiona
#     Ctrl + Enter (Windows/Linux) o Cmd + Enter (Mac).
#   • Lee los mensajes de la consola: explican cada paso y te hacen preguntas.
#   • Las figuras aparecen en el panel "Plots" y se guardan en "figuras/".
#   • Necesitas internet la primera vez (para descargar del NCBI). Después,
#     las secuencias quedan guardadas en "secuencias/" y se reutilizan.
#
###############################################################################


# =============================================================================
#  SECCIÓN 1 ✏️  CARPETA DE TRABAJO
# =============================================================================
# Escribe la ruta de la carpeta donde guardaste los dos archivos.
# Usa "/" (no "\") en la ruta. Ejemplos:
#   Windows: "C:/Users/TuNombre/Documents/actividad_secuencias"
#   Mac:     "/Users/TuNombre/Documents/actividad_secuencias"

CARPETA_TRABAJO <- "~/Desktop/Actividad_BM_DCBM/"

# Opción B (RStudio): quita el "#" de la siguiente línea para elegir la carpeta
# en una ventana.
# CARPETA_TRABAJO <- rstudioapi::selectDirectory("Elige la carpeta de la actividad")

setwd(CARPETA_TRABAJO)
source("funciones_secuencias.R", encoding = "UTF-8")
preparar_actividad()      # instala lo que falte y crea las carpetas
mostrar_menu()            # muestra los genes y especies disponibles


# =============================================================================
#  SECCIÓN 2 ✏️  PANEL DE PARÁMETROS
# =============================================================================

NOMBRE <- "Escribe aquí tu nombre"

# --- Gen (elige un número del menú) ---------------------------------------
#  1 = HBB   β-globina          5 = MB    Mioglobina
#  2 = INS   Insulina           6 = APOE  Apolipoproteína E
#  3 = CYCS  Citocromo c        7 = CCR5  Receptor CCR5
#  4 = MC1R  Receptor de melanocortina 1 (sin intrones)
GEN <- 1

# --- Especies (elige de 3 a 5 números) ------------------------------------
#  1 = Humano   2 = Chimpancé   3 = Ratón   4 = Pollo   5 = Pez cebra
#  Incluye especies cercanas (1, 2) y lejanas (4, 5) para comparar.
ESPECIES_ELEGIDAS <- c(1, 2, 3, 4, 5)

# --- Exón que vas a transcribir y traducir (módulos 3 y 4) ----------------
EXON_EJEMPLO <- 1

# --- Posición o rango para el acercamiento de la Figura 1c ----------------
#  UNA posición (se muestran 30 columnas a cada lado):
#    350      -> 350 nt después del TSS (0 = inicio del ARNm)
#    -30      -> 30 nt antes del TSS (promotor)
#    "c.20"   -> nucleótido 20 de la región codificante (la A del ATG = c.1)
#  UN RANGO (se muestra exactamente ese tramo, en filas de 100 columnas):
#    c(300, 400)        -> de +300 a +400 respecto al TSS
#    c("c.20", "c.80")  -> de c.20 a c.80 en la región codificante
#    c(-60, 20)         -> promotor y comienzo del ARNm
#  NULL = el script te lo preguntará en la consola (escribe 350, c.20 o 300:400).
COORDENADA_1C <- c(-60, 20) 

# --- Tramo para comparar ADN-ARNm-proteína (módulo 5) ---------------------
#  NULL = el script busca automáticamente el tramo con más diferencias.
#  O escribe los codones de inicio y fin (máximo 25), por ejemplo c(1, 20).
VENTANA_CODONES <- NULL

# --- Mutación que vas a rastrear (módulo 6) -------------------------------
#  Ejecuta  mostrar_mutaciones(GEN)  para ver el menú del gen elegido.
#  0 = usar tu propia mutación (definida abajo en MUTACION_PROPIA).
MUTACION_EJEMPLO <- 1

MUTACION_PROPIA <- list(
  hgvs         = "",             # opcional: notación HGVS, p. ej. "c.20A>T", "c.68_69del",
                                 # "c.30dup", "c.30_31insGGG". Si la usas, se ignoran
                                 # tipo, posicion_cds, nuevo y longitud.
  especie      = "humano",       # "humano", "chimpance", "raton", "pollo" o "pez_cebra"
  tipo         = "sustitucion",  # "sustitucion", "delecion" o "insercion"
  posicion_cds = 20,             # posición en la región codificante (A del ATG = 1)
  nuevo        = "T",            # base nueva (sustitución) o bases insertadas
  longitud     = 1               # solo para deleción: cuántas bases se eliminan
)

# --- Estructura 3D (módulo 7) ----------------------------------------------
ESPECIE_3D   <- "humano"       # especie cuya proteína verás en 3D
COMPARAR_3D  <- "raton"        # marca las diferencias con esta especie (NULL = ninguna)
COLOR_3D     <- "propiedades"  # "propiedades", "confianza" (AlphaFold) o "diferencias"

mostrar_mutaciones(GEN)


# =============================================================================
#  SECCIÓN 3  DESCARGA Y LECTURA DE LAS SECUENCIAS (no modificar)
# =============================================================================
# Descarga del NCBI la región genómica del gen (promotor + exones + intrones)
# con sus anotaciones. Si ya está en la carpeta "secuencias", la reutiliza.

datos <- cargar_datos(GEN, ESPECIES_ELEGIDAS)


# =============================================================================
#  SECCIÓN 4  MÓDULOS DE ANÁLISIS (ejecuta uno por uno y lee la consola)
# =============================================================================

# MÓDULO 1 · Alineamiento y conservación   -> Figuras 1, 1b, 2 y 3
modulo1_alineamiento(datos)

# FIGURA 1c · Acercamiento a la posición que tú elijas (repite con otras posiciones)
acercamiento(datos, COORDENADA_1C)

# MÓDULO 2 · Estructura del gen            -> Figuras 4, 5 y 5b
modulo2_estructura(datos)

# MÓDULO 3 · Transcripción de un exón      -> Figura 6
modulo3_transcripcion(datos, exon = EXON_EJEMPLO)

# MÓDULO 4 · Traducción de ese exón        -> Figura 7
modulo4_traduccion(datos, exon = EXON_EJEMPLO)

# MÓDULO 5 · Tramo variable entre especies -> Figuras 8, 8b y 8c (mapa de todos los cambios)
modulo5_comparacion(datos, ventana = VENTANA_CODONES)

# MÓDULO 6 · Efecto de una mutación        -> Figura 9
modulo6_mutacion(datos, ejemplo = MUTACION_EJEMPLO, propia = MUTACION_PROPIA)

# MÓDULO 7 · Proteína en 3D (necesita internet) -> visor rotable + Figura 11
#  Se abre en el panel "Viewer" de RStudio: gírala con el mouse.
#  También queda en la carpeta "estructuras/" (HTML, .pdb y scripts de PyMOL/ChimeraX).
modulo7_estructura3d(datos, especie = ESPECIE_3D, comparar_con = COMPARAR_3D, color = COLOR_3D)


# =============================================================================
#  SECCIÓN 4b  HERRAMIENTAS DE BÚSQUEDA (úsalas para responder las preguntas 🔎)
# =============================================================================
# Cambia los valores y ejecuta las que necesites:
# buscar_secuencia(datos, "GTAAGT")          # busca una secuencia (admite IUPAC: R, Y, W, N...)
# buscar_secuencia(datos, "TATAWAW", max_dif = 1)   # permite 1 diferencia
# ver_region(datos, -60, 20)                 # muestra un tramo (posiciones respecto al TSS)
# uso_codones(datos)                         # contenido GC por posición del codón
# perfil_gc(datos)                           # contenido GC e islas CpG a lo largo del gen
# acercamiento(datos, "c.100")               # otro acercamiento (Figura 1c)
# alinear_exon(datos, exon = 2)              # un exón con sus bordes intrónicos en todas las especies


# =============================================================================
#  SECCIÓN 5 ✏️  MÓDULO DE RESPUESTAS
# =============================================================================
# Escribe las secuencias entre comillas (en ADN: A, C, G, T).
# Usa las figuras y las secuencias anotadas (carpeta resultados/) como ayuda.
# Si tu respuesta es correcta, el script la buscará en todas las especies y
# mostrará dónde aparece (Figura 10).

RESPUESTAS <- list(
  codon_inicio  = "",            # codón donde empieza la traducción
  codones_stop  = c("", "", ""), # los tres codones de parada
  caja_TATA     = "",            # secuencia consenso de la caja TATA del promotor
  caja_CCAAT    = "",            # otra caja del promotor, 60-100 nt antes del TSS
  sitio_donador = "",            # dinucleótido con el que EMPIEZA un intrón
  sitio_aceptor = "",            # dinucleótido con el que TERMINA un intrón
  senal_poliA   = ""             # señal de poliadenilación cerca del final del ARNm
)

revisar_respuestas(datos, RESPUESTAS)


# =============================================================================
#  SECCIÓN 6  INFORME FINAL
# =============================================================================
# Genera un resumen en resultados/ con tablas, mutación, respuestas y figuras.

generar_informe(datos, NOMBRE)


###############################################################################
#  PARA IR MÁS ALLÁ (opcional)
#  • Cambia GEN y repite: ¿los genes muy conservados (CYCS) se comportan igual
#    que los que cambian rápido?
#  • Compara HBB (3 exones) con MC1R (sin intrones) en la Figura 4.
#  • Con APOE, busca el codón 130 en la Figura 8 (VENTANA_CODONES <- c(120, 140)):
#    ¿qué aminoácido tienen el chimpancé y el ratón en la posición del alelo ε4?
#
#  PARA EL DOCENTE
#  • Para trabajar sin internet en clase, ejecuta una vez (con internet):
#        source("funciones_secuencias.R"); preparar_actividad(); descargar_todo()
#    y comparte la carpeta "secuencias" con los estudiantes.
###############################################################################
