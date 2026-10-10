# =============================================================================
#  funciones_secuencias.R
#  Motor de la actividad "Del gen a la proteína en cinco especies"
#  Genética Molecular Evolutiva Humana · Universidad del Tolima
#
#  >>> NO ES NECESARIO MODIFICAR ESTE ARCHIVO <<<
#  El estudiante trabaja en "actividad_secuencias.R", que carga este archivo.
#
#  Contenido
#    0. Utilidades y paquetes
#    1. Catálogo de genes, especies y mutaciones
#    2. Descarga desde NCBI (con copia local en la carpeta "secuencias")
#    3. Lectura de archivos GenBank y construcción del modelo del gen
#    4. Alineamientos e identidad
#    5. Módulos de la actividad (1 a 6) y sus figuras
#    6. Revisión de respuestas
#    7. Informe final
# =============================================================================


# -----------------------------------------------------------------------------
# 0. UTILIDADES Y PAQUETES
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

PAQUETES_CRAN <- c("rentrez", "ggplot2", "patchwork")
PAQUETES_BIOC <- c("Biostrings", "msa")

# Instala lo que falte, carga los paquetes y crea las carpetas de trabajo.
preparar_actividad <- function(carpeta = getwd()) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  faltan <- PAQUETES_CRAN[!vapply(PAQUETES_CRAN, requireNamespace,
                                  logical(1), quietly = TRUE)]
  if (length(faltan)) {
    message("Instalando paquetes de CRAN: ", paste(faltan, collapse = ", "))
    install.packages(faltan)
  }
  faltan_b <- PAQUETES_BIOC[!vapply(PAQUETES_BIOC, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(faltan_b)) {
    message("Instalando paquetes de Bioconductor: ",
            paste(faltan_b, collapse = ", "), " (puede tardar unos minutos)")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(faltan_b, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages({
    library(Biostrings)
    library(ggplot2)
    library(patchwork)
  })
  for (d in c("secuencias", "figuras", "resultados"))
    dir.create(file.path(carpeta, d), showWarnings = FALSE)
  .estado$carpeta <- normalizePath(carpeta)
  cat("\n✔ Paquetes cargados y carpetas listas en:\n  ", .estado$carpeta, "\n",
      "   secuencias/  -> archivos descargados del NCBI\n",
      "   figuras/     -> gráficos generados (PNG)\n",
      "   resultados/  -> tablas e informe final\n\n", sep = "")
  invisible(TRUE)
}

.estado <- new.env()
.estado$carpeta <- getwd()

ruta <- function(...) file.path(.estado$carpeta, ...)

caja <- function(titulo, ...) {
  linea <- strrep("─", 72)
  cat("\n", linea, "\n  ", titulo, "\n", linea, "\n", sep = "")
  txt <- c(...)
  if (length(txt)) cat(paste0("  ", txt), sep = "\n")
  cat("\n")
}

preguntas <- function(...) {
  cat("  ✎ Preguntas para tu informe:\n")
  q <- c(...)
  cat(paste0("    ", seq_along(q), ". ", q), sep = "\n")
  cat("\n")
}

letras <- function(s) strsplit(s, "")[[1]]

revcomp <- function(s) as.character(reverseComplement(DNAString(s)))

traducir <- function(s) {
  n3 <- floor(nchar(s) / 3) * 3
  if (n3 < 3) return("")
  as.character(suppressWarnings(
    translate(DNAString(substr(s, 1, n3)), no.init.codon = TRUE, if.fuzzy.codon = "X")))
}

# Traduce hasta el primer codón de parada (incluido como "*").
traducir_hasta_stop <- function(s) {
  p <- traducir(s)
  k <- regexpr("\\*", p)
  if (k > 0) substr(p, 1, k) else p
}

aa3 <- function(a) {
  a <- toupper(a)
  out <- unname(Biostrings::AMINO_ACID_CODE[a])
  out[a == "*"] <- "Ter"
  out[a == "-"] <- "---"
  out[is.na(out)] <- "Xaa"
  out
}

guardar_figura <- function(p, nombre, ancho, alto) {
  archivo <- ruta("figuras", paste0(.estado$gen, "_", nombre, ".png"))
  ggsave(archivo, p, width = ancho, height = alto, dpi = 150, bg = "white")
  .estado$figuras <- unique(c(.estado$figuras, basename(archivo)))
  cat("  \U0001F5BC  Figura guardada: figuras/", basename(archivo), "\n", sep = "")
  invisible(archivo)
}

mostrar <- function(p) {
  print(p)
  invisible(p)
}

# Colores fijos para cada nucleótido (los mismos en todas las figuras)
COL_NT <- c(A = "#2E9E44", C = "#1F77B4", G = "#F2A900", T = "#D62728",
            U = "#D62728", "-" = "#E9E9E9", N = "#9E9E9E")

COL_REG <- c("promotor" = "#CFC6B0", "5′UTR" = "#8DB8D0",
             "CDS" = "#1F4E79", "intrón" = "#F1ECE2", "3′UTR" = "#8DB8D0",
             "intergénico" = "#FFFFFF")

tema_base <- function(tam = 12) {
  theme_minimal(base_size = tam) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 1. CATÁLOGO DE GENES, ESPECIES Y MUTACIONES
# -----------------------------------------------------------------------------

ESPECIES_INFO <- data.frame(
  clave      = c("humano", "chimpance", "raton", "pollo", "pez_cebra"),
  nombre     = c("Humano", "Chimpancé", "Ratón", "Pollo", "Pez cebra"),
  cientifico = c("Homo sapiens", "Pan troglodytes", "Mus musculus",
                 "Gallus gallus", "Danio rerio"),
  grupo      = c("Primate", "Primate", "Roedor", "Ave", "Pez óseo"),
  # Tiempo aproximado desde el ancestro común con el humano (millones de años)
  divergencia_ma = c(0, 6.5, 90, 320, 435),
  taxon      = c(9606, 9598, 10090, 9031, 7955),   # identificador taxonómico del NCBI
  stringsAsFactors = FALSE
)

CATALOGO <- list(
  list(
    clave = "HBB", nombre = "β-globina (hemoglobina, cadena beta)",
    simbolos = list(humano = "HBB", chimpance = "HBB",
                    raton = c("Hbb-bs", "Hbb-b1", "Hbb-bt"),
                    pollo = c("HBB", "HBBA"),
                    pez_cebra = c("hbba1", "hbba1.1", "hbba2")),
    resumen = c("Gen corto con 3 exones y 2 intrones. Transporta oxígeno en la sangre.",
                "Sus mutaciones causan anemia falciforme y talasemias; el alelo HbS",
                "se mantiene alto en regiones con malaria (ventaja del heterocigoto)."),
    mutaciones = list(
      list(etiqueta = "Anemia falciforme (HbS): c.20A>T, p.Glu7Val (Glu6Val en numeración clásica)",
           tipo = "sustitucion", posicion = 20, ref = "A", nuevo = "T"),
      list(etiqueta = "Hemoglobina C: c.19G>A, p.Glu7Lys (Glu6Lys clásica)",
           tipo = "sustitucion", posicion = 19, ref = "G", nuevo = "A"),
      list(etiqueta = "β0-talasemia, codón 39 sin sentido: c.118C>T, p.Gln40Ter",
           tipo = "sustitucion", posicion = 118, ref = "C", nuevo = "T"))),
  list(
    clave = "INS", nombre = "Insulina",
    simbolos = list(humano = "INS", chimpance = "INS", raton = c("Ins2", "Ins1"),
                    pollo = "INS", pez_cebra = c("ins", "insa")),
    resumen = c("Hormona que regula la glucosa en sangre; gen de 3 exones.",
                "El ratón tiene dos copias (parálogos): Ins2 conserva los dos intrones",
                "e Ins1 perdió uno. Mutaciones en INS causan diabetes neonatal."),
    mutaciones = list()),
  list(
    clave = "CYCS", nombre = "Citocromo c",
    simbolos = list(humano = "CYCS", chimpance = "CYCS", raton = "Cycs",
                    pollo = "CYCS", pez_cebra = c("cycsb", "cycs", "cycsa")),
    resumen = c("Proteína de la cadena respiratoria y de la apoptosis, muy conservada.",
                "Fitch y Margoliash (1967) construyeron uno de los primeros árboles",
                "filogenéticos moleculares con su secuencia de aminoácidos."),
    mutaciones = list()),
  list(
    clave = "MC1R", nombre = "Receptor de melanocortina 1",
    simbolos = list(humano = "MC1R", chimpance = "MC1R", raton = "Mc1r",
                    pollo = "MC1R", pez_cebra = "mc1r"),
    resumen = c("Controla el tipo de melanina (eumelanina o feomelanina).",
                "Su región codificante NO tiene intrones: sirve de contraste.",
                "Variantes como R151C se asocian con pelo rojo, piel clara y melanoma."),
    mutaciones = list(
      list(etiqueta = "Pelo rojo / riesgo de melanoma: c.451C>T, p.Arg151Cys (R151C)",
           tipo = "sustitucion", posicion = 451, ref = "C", nuevo = "T"),
      list(etiqueta = "Pelo rojo: c.478C>T, p.Arg160Trp (R160W)",
           tipo = "sustitucion", posicion = 478, ref = "C", nuevo = "T"))),
  list(
    clave = "MB", nombre = "Mioglobina",
    simbolos = list(humano = "MB", chimpance = "MB", raton = "Mb",
                    pollo = "MB", pez_cebra = "mb"),
    resumen = c("Almacena oxígeno en el músculo. Comparte un ancestro lejano con la",
                "β-globina (son parálogos): compara su estructura de 3 exones con HBB."),
    mutaciones = list()),
  list(
    clave = "APOE", nombre = "Apolipoproteína E",
    simbolos = list(humano = "APOE", chimpance = "APOE", raton = "Apoe",
                    pollo = "APOE", pez_cebra = c("apoea", "apoeb")),
    resumen = c("Transporta colesterol y lípidos. Los alelos ε2, ε3 y ε4 difieren en los",
                "codones 130 y 176 (112 y 158 en numeración clásica); ε4 aumenta el riesgo",
                "de Alzheimer. ε4 es el alelo ancestral: otros primates tienen Arg en 112.",
                "Si una especie no aparece, puede no tener el gen anotado (¡o no tenerlo!)."),
    mutaciones = list(
      list(etiqueta = "Alelo ε4 (desde ε3): c.388T>C, p.Cys130Arg (Cys112Arg clásica)",
           tipo = "sustitucion", posicion = 388, ref = "T", nuevo = "C"),
      list(etiqueta = "Alelo ε2 (desde ε3): c.526C>T, p.Arg176Cys (Arg158Cys clásica)",
           tipo = "sustitucion", posicion = 526, ref = "C", nuevo = "T"))),
  list(
    clave = "CCR5", nombre = "Receptor de quimiocinas CCR5",
    simbolos = list(humano = "CCR5", chimpance = "CCR5", raton = "Ccr5",
                    pollo = "CCR5", pez_cebra = "ccr5"),
    resumen = c("Correceptor que usa el VIH para entrar a los linfocitos.",
                "La deleción Δ32 (32 pb) produce una proteína truncada: los homocigotos",
                "son resistentes al VIH. Su origen y frecuencia en Europa se debaten.",
                "Si una especie no aparece, puede no tener un ortólogo claro anotado."),
    mutaciones = list(
      list(etiqueta = "CCR5-Δ32: c.554_585del (deleción de 32 pb, cambio de marco)",
           tipo = "delecion", posicion = 554, longitud = 32)))
)
names(CATALOGO) <- vapply(CATALOGO, `[[`, "", "clave")

MUTACIONES_AUTO <- list(
  list(etiqueta = "Automática: mutación SILENCIOSA (cambia la 3.ª base de un codón)",
       tipo = "auto_silenciosa"),
  list(etiqueta = "Automática: mutación SIN SENTIDO (crea un codón de parada prematuro)",
       tipo = "auto_sin_sentido"),
  list(etiqueta = "Automática: CAMBIO DE MARCO (deleción de 1 nucleótido en el codón 10)",
       tipo = "auto_cambio_marco")
)

menu_mutaciones <- function(gen) c(CATALOGO[[gen]]$mutaciones, MUTACIONES_AUTO)

resolver_gen <- function(gen) {
  if (is.numeric(gen)) {
    if (gen < 1 || gen > length(CATALOGO)) stop("GEN debe ser un número entre 1 y ", length(CATALOGO))
    return(names(CATALOGO)[gen])
  }
  g <- toupper(gen)
  if (!g %in% names(CATALOGO)) stop("El gen '", gen, "' no está en el menú.")
  g
}

resolver_especies <- function(esp) {
  if (is.numeric(esp)) {
    if (any(esp < 1 | esp > nrow(ESPECIES_INFO))) stop("Especies: usa números del 1 al ", nrow(ESPECIES_INFO))
    return(ESPECIES_INFO$clave[sort(unique(esp))])
  }
  esp <- tolower(esp)
  malas <- setdiff(esp, ESPECIES_INFO$clave)
  if (length(malas)) stop("Especies no reconocidas: ", paste(malas, collapse = ", "))
  ESPECIES_INFO$clave[ESPECIES_INFO$clave %in% esp]
}

nombre_esp <- function(clave) ESPECIES_INFO$nombre[match(clave, ESPECIES_INFO$clave)]

mostrar_menu <- function() {
  caja("MENÚ DE GENES (usa el número en GEN)")
  for (i in seq_along(CATALOGO)) {
    g <- CATALOGO[[i]]
    cat(sprintf("  %d = %-5s %s\n", i, g$clave, g$nombre))
    cat(paste0("            ", g$resumen), sep = "\n")
  }
  caja("MENÚ DE ESPECIES (usa los números en ESPECIES_ELEGIDAS)")
  for (i in seq_len(nrow(ESPECIES_INFO)))
    cat(sprintf("  %d = %-10s (%s) - %s, ~%s millones de años desde el ancestro común con el humano\n",
                i, ESPECIES_INFO$nombre[i], ESPECIES_INFO$cientifico[i],
                ESPECIES_INFO$grupo[i], ESPECIES_INFO$divergencia_ma[i]))
  cat("\n  Para ver las mutaciones de un gen: mostrar_mutaciones(GEN)\n\n")
  invisible(NULL)
}

mostrar_mutaciones <- function(gen) {
  gen <- resolver_gen(gen)
  m <- menu_mutaciones(gen)
  caja(paste0("MENÚ DE MUTACIONES PARA ", gen, " (usa el número en MUTACION_EJEMPLO)"))
  cat("  0 = Mutación propia (defínela en MUTACION_PROPIA)\n")
  for (i in seq_along(m)) cat(sprintf("  %d = %s\n", i, m[[i]]$etiqueta))
  cat("\n")
  invisible(m)
}


# -----------------------------------------------------------------------------
# 2. DESCARGA DESDE NCBI
#    Se descarga la REGIÓN GENÓMICA del gen (con su promotor), en formato GenBank,
#    que incluye la secuencia y las anotaciones (ARNm, CDS = región codificante).
# -----------------------------------------------------------------------------

FLANCO_5 <- 1000   # nt antes del gen (promotor)
FLANCO_3 <- 500    # nt después del gen

archivo_gb   <- function(gen, esp) ruta("secuencias", paste0(gen, "_", esp, ".gb"))
archivo_meta <- function(gen, esp) ruta("secuencias", paste0(gen, "_", esp, "_info.txt"))

buscar_gen_ncbi <- function(simbolos, cientifico) {
  for (s in simbolos) {
    termino <- sprintf('%s[sym] AND "%s"[orgn] AND alive[prop]', s, cientifico)
    r <- rentrez::entrez_search(db = "gene", term = termino, retmax = 5)
    Sys.sleep(0.4)
    if (length(r$ids) == 0) next
    res <- rentrez::entrez_summary(db = "gene", id = r$ids)
    Sys.sleep(0.4)
    if (!inherits(res, "esummary_list")) res <- list(res)
    nombres <- vapply(res, function(x) x$name %||% "", "")
    k <- which(tolower(nombres) == tolower(s))[1]
    if (is.na(k)) k <- 1
    return(list(id = r$ids[k], resumen = res[[k]], simbolo = nombres[k]))
  }
  NULL
}

descargar_gen <- function(gen, esp) {
  info <- ESPECIES_INFO[ESPECIES_INFO$clave == esp, ]
  simbolos <- CATALOGO[[gen]]$simbolos[[esp]] %||% CATALOGO[[gen]]$simbolos$humano
  cat(sprintf("  Buscando %s en %s (%s) ... ", gen, info$nombre, info$cientifico))
  hallado <- buscar_gen_ncbi(simbolos, info$cientifico)
  if (is.null(hallado)) {
    cat("no encontrado.\n")
    return(FALSE)
  }
  gi <- hallado$resumen$genomicinfo
  if (is.data.frame(gi)) gi <- gi[1, , drop = FALSE] else if (is.list(gi) && length(gi)) gi <- gi[[1]]
  if (is.null(gi$chraccver)) {
    cat("sin coordenadas genómicas en NCBI.\n")
    return(FALSE)
  }
  acc <- gi$chraccver
  ini <- as.numeric(gi$chrstart) + 1   # NCBI usa coordenadas desde 0
  fin <- as.numeric(gi$chrstop) + 1
  hebra <- if (ini > fin) "-" else "+"
  a <- min(ini, fin); b <- max(ini, fin)
  if (hebra == "+") { a <- a - FLANCO_5; b <- b + FLANCO_3 } else { a <- a - FLANCO_3; b <- b + FLANCO_5 }
  a <- max(1, a)
  txt <- rentrez::entrez_fetch(db = "nuccore", id = acc, rettype = "gbwithparts",
                               retmode = "text", seq_start = a, seq_stop = b, strand = 1)
  Sys.sleep(0.4)
  if (!grepl("^LOCUS", txt)) {
    cat("la descarga no devolvió un registro GenBank.\n")
    return(FALSE)
  }
  writeLines(txt, archivo_gb(gen, esp))
  meta <- data.frame(gen = gen, especie = esp, simbolo = hallado$simbolo,
                     gene_id = hallado$id, descripcion = hallado$resumen$description %||% "",
                     accesion = acc, inicio = a, fin = b, hebra = hebra,
                     fecha = format(Sys.Date()), stringsAsFactors = FALSE)
  write.dcf(meta, archivo_meta(gen, esp))
  cat(sprintf("listo (%s, Gene ID %s, %s:%s-%s, hebra %s)\n",
              hallado$simbolo, hallado$id, acc, a, b, hebra))
  TRUE
}

# Para el docente: descarga todos los genes y especies una sola vez y
# comparte la carpeta "secuencias" con los estudiantes (trabajo sin internet).
descargar_todo <- function(genes = names(CATALOGO), especies = ESPECIES_INFO$clave) {
  for (g in genes) for (e in especies) {
    if (file.exists(archivo_gb(g, e))) next
    ok <- tryCatch(descargar_gen(g, e), error = function(err) {
      cat("error:", conditionMessage(err), "\n"); FALSE })
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# 3. LECTURA DE GENBANK Y MODELO DEL GEN
# -----------------------------------------------------------------------------

leer_genbank <- function(texto) {
  if (length(texto) > 1) texto <- paste(texto, collapse = "\n")
  lin <- strsplit(texto, "\r?\n")[[1]]
  fin_reg <- which(lin == "//")
  if (length(fin_reg)) lin <- lin[seq_len(fin_reg[1] - 1)]
  campo <- function(clave) {
    i <- grep(paste0("^", clave, "\\s"), lin)[1]
    if (is.na(i)) "" else trimws(sub(paste0("^", clave), "", lin[i]))
  }
  i_f <- grep("^FEATURES", lin)[1]
  i_o <- grep("^ORIGIN", lin)[1]
  if (is.na(i_f) || is.na(i_o)) stop("El archivo GenBank no tiene anotaciones o secuencia.")
  sig <- grep("^\\S", lin)
  fin_f <- sig[sig > i_f][1] - 1
  fl <- lin[(i_f + 1):fin_f]
  ini_feat <- grep("^ {5}\\S", fl)
  features <- vector("list", length(ini_feat))
  for (k in seq_along(ini_feat)) {
    a <- ini_feat[k]
    b <- if (k < length(ini_feat)) ini_feat[k + 1] - 1 else length(fl)
    bloque <- fl[a:b]
    tipo <- trimws(substr(bloque[1], 6, 21))
    cuerpo <- trimws(substring(bloque, 22))
    es_q <- startsWith(cuerpo, "/")
    pq <- which(es_q)[1]
    if (is.na(pq)) pq <- length(cuerpo) + 1
    ubic <- paste(cuerpo[seq_len(pq - 1)], collapse = "")
    qual <- list()
    if (pq <= length(cuerpo)) {
      resto <- cuerpo[pq:length(cuerpo)]
      grupos <- split(resto, cumsum(startsWith(resto, "/")))
      for (p in grupos) {
        m <- regmatches(p[1], regexec("^/([^=]+)(=(.*))?$", p[1]))[[1]]
        if (!length(m)) next
        nom <- m[2]
        val <- m[4]
        if (length(p) > 1)
          val <- paste(c(val, p[-1]), collapse = if (nom == "translation") "" else " ")
        val <- gsub('^"|"$', "", val)
        qual[[nom]] <- c(qual[[nom]], val)
      }
    }
    features[[k]] <- list(tipo = tipo, ubicacion = ubic, qual = qual)
  }
  sec <- toupper(gsub("[^A-Za-z]", "", paste(lin[(i_o + 1):length(lin)], collapse = "")))
  list(locus = campo("LOCUS"), definicion = campo("DEFINITION"),
       accesion = campo("ACCESSION"), version = campo("VERSION"),
       features = features, secuencia = sec)
}

parsear_ubicacion <- function(loc) {
  loc <- gsub("\\s", "", loc)
  loc <- gsub("one-of\\(([0-9]+)(,[0-9]+)*\\)", "\\1", loc)
  remota <- grepl(":", loc)
  hebra <- if (grepl("complement\\(", loc)) "-" else "+"
  partes <- strsplit(gsub("complement\\(|join\\(|order\\(|\\)", "", loc), ",")[[1]]
  partes <- partes[!grepl(":", partes) & nzchar(partes)]
  if (!length(partes))
    return(list(rangos = data.frame(inicio = integer(), fin = integer()),
                hebra = hebra, remota = TRUE))
  ini <- fin <- integer(length(partes))
  for (i in seq_along(partes)) {
    n <- as.integer(regmatches(partes[i], gregexpr("[0-9]+", partes[i]))[[1]])
    if (grepl("\\.\\.", partes[i]) && length(n) >= 2) { ini[i] <- n[1]; fin[i] <- n[2] }
    else { ini[i] <- n[1]; fin[i] <- n[1] }
  }
  r <- data.frame(inicio = pmin(ini, fin), fin = pmax(ini, fin))
  r <- r[order(r$inicio), , drop = FALSE]
  rownames(r) <- NULL
  list(rangos = r, hebra = hebra, remota = remota,
       parcial = grepl("[<>]", loc))
}

posiciones <- function(r) unlist(mapply(seq, r$inicio, r$fin, SIMPLIFY = FALSE))

# Construye el modelo del gen: secuencia orientada 5'->3', exones, intrones,
# CDS, sitio de inicio de la transcripción (TSS), codón de inicio, parada, etc.
construir_modelo <- function(gb, meta, promotor = 200) {
  L <- nchar(gb$secuencia)
  feats <- gb$features
  tipos <- vapply(feats, `[[`, "", "tipo")
  gen_q <- vapply(feats, function(f) tolower(f$qual$gene[1] %||% ""), "")
  ub <- lapply(feats, function(f) parsear_ubicacion(f$ubicacion))

  # Seguridad: si las coordenadas vienen en el sistema del cromosoma, se ajustan.
  mx <- max(c(0, unlist(lapply(ub, function(u) u$rangos$fin))))
  if (mx > L && !is.null(meta$inicio) && as.numeric(meta$inicio) > 1) {
    off <- as.numeric(meta$inicio) - 1
    ub <- lapply(ub, function(u) { u$rangos$inicio <- u$rangos$inicio - off
                                   u$rangos$fin <- u$rangos$fin - off; u })
  }
  dentro <- vapply(ub, function(u) nrow(u$rangos) > 0 && !u$remota &&
                     min(u$rangos$inicio) >= 1 && max(u$rangos$fin) <= L, TRUE)

  sim <- tolower(meta$simbolo)
  candidatos <- function(tipo) {
    idx <- which(tipos == tipo & gen_q == sim & dentro)
    if (!length(idx)) {
      # Respaldo: elementos que cubren el centro de la región descargada
      centro <- round(L / 2)
      idx <- which(tipos == tipo & dentro & vapply(ub, function(u)
        min(u$rangos$inicio) <= centro + L / 4 && max(u$rangos$fin) >= centro - L / 4, TRUE))
    }
    idx
  }
  # Preferencia: MANE Select (transcrito canónico) > NM_/NP_ > XM_/XP_ > más largo
  rango_pref <- function(idx, campo, prefs) {
    if (!length(idx)) return(idx)
    ids <- vapply(idx, function(i) feats[[i]]$qual[[campo]][1] %||% "", "")
    mane <- vapply(idx, function(i) any(grepl("MANE Select", feats[[i]]$qual$tag %||% "")), TRUE)
    largo <- vapply(idx, function(i) sum(ub[[i]]$rangos$fin - ub[[i]]$rangos$inicio + 1), 1)
    idx[order(!mane, match(substr(ids, 1, 3), prefs, nomatch = 99), -largo)]
  }
  im <- rango_pref(candidatos("mRNA"), "transcript_id", c("NM_", "XM_"))
  ic <- rango_pref(candidatos("CDS"), "protein_id", c("NP_", "XP_"))
  if (!length(ic)) stop("No se encontró una región codificante (CDS) anotada para ", meta$simbolo)

  elegido_m <- NA; elegido_c <- NA
  for (i in im) {
    pm <- posiciones(ub[[i]]$rangos)
    for (j in ic) {
      if (ub[[j]]$hebra != ub[[i]]$hebra) next
      if (all(posiciones(ub[[j]]$rangos) %in% pm)) { elegido_m <- i; elegido_c <- j; break }
    }
    if (!is.na(elegido_m)) break
  }
  if (is.na(elegido_c)) elegido_c <- ic[1]
  cds <- ub[[elegido_c]]$rangos
  hebra <- ub[[elegido_c]]$hebra
  sin_utr <- is.na(elegido_m)
  exones <- if (sin_utr) cds else ub[[elegido_m]]$rangos
  qc <- feats[[elegido_c]]$qual

  sec <- gb$secuencia
  if (hebra == "-") {
    sec <- revcomp(sec)
    voltear <- function(r) {
      out <- data.frame(inicio = L - r$fin + 1, fin = L - r$inicio + 1)
      out[rev(seq_len(nrow(out))), , drop = FALSE]
    }
    cds <- voltear(cds); exones <- voltear(exones)
  }
  rownames(cds) <- NULL; rownames(exones) <- NULL
  nt <- letras(sec)
  mapa_cds <- posiciones(cds)
  mapa_arn <- posiciones(exones)
  cs <- as.integer(qc$codon_start[1] %||% 1)
  cds_sec <- paste(nt[mapa_cds], collapse = "")
  cds_util <- substring(cds_sec, cs)
  prot <- traducir(cds_util)
  prot_ss <- sub("\\*$", "", prot)
  tr_ncbi <- gsub("\\s", "", qc$translation[1] %||% "")
  valido <- nzchar(tr_ncbi) &&
    (prot_ss == tr_ncbi || substring(prot_ss, 2) == substring(tr_ncbi, 2))

  tss <- exones$inicio[1]; tes <- exones$fin[nrow(exones)]
  atg <- cds$inicio[1]; fin_cds <- cds$fin[nrow(cds)]
  intrones <- if (nrow(exones) > 1)
    data.frame(numero = seq_len(nrow(exones) - 1),
               inicio = exones$fin[-nrow(exones)] + 1, fin = exones$inicio[-1] - 1)
  else data.frame(numero = integer(), inicio = integer(), fin = integer())
  exones$numero <- seq_len(nrow(exones))

  region <- rep("intergénico", L)
  if (tss > 1) region[max(1, tss - promotor):(tss - 1)] <- "promotor"
  region[tss:tes] <- "intrón"
  region[mapa_arn] <- ifelse(mapa_arn < atg, "5′UTR", "3′UTR")
  region[mapa_cds] <- "CDS"

  # Señal de poliadenilación cerca del final del transcrito
  zona <- substr(sec, max(1, tes - 45), tes)
  hits <- gregexpr("AATAAA|ATTAAA", zona)[[1]]
  poliA <- if (hits[1] > 0) max(1, tes - 45) + tail(hits, 1) - 1 else NA

  list(gen = meta$gen, especie = meta$especie, simbolo = meta$simbolo,
       gene_id = meta$gene_id, accesion = meta$accesion,
       descripcion = meta$descripcion,
       transcrito = if (sin_utr) "(sin ARNm anotado)" else feats[[elegido_m]]$qual$transcript_id[1] %||% "",
       proteina_id = qc$protein_id[1] %||% "",
       hebra_original = hebra, sec = sec, nt = nt, L = L,
       exones = exones, intrones = intrones, cds = cds,
       mapa_cds = mapa_cds, mapa_arn = mapa_arn,
       tss = tss, tes = tes, atg = atg, fin_cds = fin_cds, poliA = poliA,
       codon_inicio = cs,
       cds_sec = cds_util, arnm = paste(nt[mapa_arn], collapse = ""),
       proteina = prot_ss, prot_ncbi = tr_ncbi, valido = valido,
       sin_utr = sin_utr, region = region)
}

# Carga (o descarga) los genes elegidos y construye sus modelos.
cargar_datos <- function(gen, especies, forzar_descarga = FALSE) {
  gen <- resolver_gen(gen)
  especies <- resolver_especies(especies)
  .estado$gen <- gen
  .estado$figuras <- character()
  caja(paste0("Cargando el gen ", gen, " - ", CATALOGO[[gen]]$nombre),
       CATALOGO[[gen]]$resumen)
  modelos <- list()
  for (e in especies) {
    if (forzar_descarga || !file.exists(archivo_gb(gen, e))) {
      ok <- tryCatch(descargar_gen(gen, e), error = function(err) {
        cat("\n  ✘ No se pudo conectar con NCBI:", conditionMessage(err), "\n",
            "    Revisa tu conexión a internet o pide al docente la carpeta 'secuencias'.\n")
        FALSE })
      if (!ok) next
    }
    meta <- as.list(read.dcf(archivo_meta(gen, e))[1, ])
    gb <- leer_genbank(readLines(archivo_gb(gen, e), warn = FALSE))
    m <- tryCatch(construir_modelo(gb, meta), error = function(err) {
      cat("  ✘", nombre_esp(e), ":", conditionMessage(err), "\n"); NULL })
    if (!is.null(m)) modelos[[e]] <- m
  }
  if (length(modelos) < 2) stop("Se necesitan al menos 2 especies con el gen disponible.")
  if (length(modelos) < 3) message("Aviso: solo hay ", length(modelos), " especies; la actividad pide al menos 3.")

  datos <- new.env()
  datos$gen <- gen
  datos$modelos <- modelos
  datos$especies <- names(modelos)
  datos$ref <- if ("humano" %in% names(modelos)) "humano" else names(modelos)[1]

  tabla <- data.frame(
    Especie = nombre_esp(names(modelos)),
    Simbolo = vapply(modelos, `[[`, "", "simbolo"),
    Accesion = vapply(modelos, `[[`, "", "accesion"),
    Transcrito = vapply(modelos, `[[`, "", "transcrito"),
    `Gen (nt)` = vapply(modelos, function(m) m$tes - m$tss + 1, 1),
    Exones = vapply(modelos, function(m) nrow(m$exones), 1),
    `CDS (nt)` = vapply(modelos, function(m) nchar(m$cds_sec), 1),
    Aminoacidos = vapply(modelos, function(m) nchar(m$proteina), 1),
    `Validación` = ifelse(vapply(modelos, `[[`, TRUE, "valido"), "✔", "revisar"),
    check.names = FALSE, row.names = NULL)
  datos$tabla <- tabla
  cat("\n")
  print(tabla, row.names = FALSE)
  cat("\n  Validación ✔ = la proteína traducida por R coincide con la del NCBI.\n")
  cat("  Especie de referencia para las comparaciones:", nombre_esp(datos$ref), "\n\n")
  write.csv(tabla, ruta("resultados", paste0(gen, "_tabla_especies.csv")), row.names = FALSE)
  datos
}


# -----------------------------------------------------------------------------
# 4. ALINEAMIENTOS E IDENTIDAD
# -----------------------------------------------------------------------------

alinear <- function(seqs) {
  if (requireNamespace("msa", quietly = TRUE)) {
    utils::capture.output(a <- suppressMessages(suppressWarnings(
      msa::msa(seqs, method = "ClustalOmega", order = "input"))))
    out <- as.character(a)
  } else if (requireNamespace("DECIPHER", quietly = TRUE)) {
    out <- as.character(DECIPHER::AlignSeqs(seqs, verbose = FALSE))
  } else stop("Instala el paquete 'msa': BiocManager::install('msa')")
  out[names(seqs)]
}

# Alineamiento de la región codificante guiado por la proteína (codón a codón),
# para que los huecos nunca rompan el marco de lectura.
alinear_codones <- function(cds, aln_prot) {
  out <- character(length(cds))
  for (i in seq_along(cds)) {
    s <- cds[i]
    n <- floor(nchar(s) / 3)
    codones <- substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
    aa <- letras(aln_prot[i])
    k <- 0
    res <- character(length(aa))
    for (j in seq_along(aa)) {
      if (aa[j] == "-") res[j] <- "---" else { k <- k + 1; res[j] <- codones[k] %||% "NNN" }
    }
    out[i] <- paste(res, collapse = "")
  }
  names(out) <- names(cds)
  out
}

matriz_aln <- function(aln) do.call(rbind, lapply(aln, letras))

identidad <- function(aln) {
  M <- matriz_aln(aln)
  n <- nrow(M)
  I <- matrix(100, n, n, dimnames = list(rownames(M), rownames(M)))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    cols <- M[i, ] != "-" | M[j, ] != "-"
    I[i, j] <- 100 * mean(M[i, cols] == M[j, cols])
  }
  I
}

# Posición (en la secuencia sin huecos) que ocupa cada columna del alineamiento
columnas_a_posicion <- function(fila) {
  p <- cumsum(fila != "-")
  p[fila == "-"] <- NA
  p
}


# -----------------------------------------------------------------------------
# 5. MÓDULOS DE LA ACTIVIDAD
# -----------------------------------------------------------------------------

orden_especies <- function(datos) {
  e <- datos$especies
  c(datos$ref, setdiff(e[order(match(e, ESPECIES_INFO$clave))], datos$ref))
}

# ---------- MÓDULO 1: alineamiento comparativo -------------------------------
modulo1_alineamiento <- function(datos, promotor = 200, ventana = 30, max_gen = 25000) {
  caja("MÓDULO 1 · Alineamiento del mismo gen en varias especies",
       "Se alinean tres versiones del gen: (1) el gen completo con su promotor,",
       "(2) la región codificante (CDS) y (3) la proteína.",
       "Cada nucleótido tiene SIEMPRE el mismo color: A verde, C azul, G amarillo, T rojo.",
       "Las bandas verticales de un mismo color son regiones conservadas.")
  esp <- orden_especies(datos)
  mods <- datos$modelos[esp]
  ref <- datos$ref
  mr <- mods[[ref]]
  largo_max <- max(vapply(mods, function(m) m$tes - m$tss + 1, 1))
  modo <- if (largo_max > max_gen) "cds" else "gen"

  cat("  Alineando proteínas y regiones codificantes...\n")
  aln_prot <- alinear(AAStringSet(vapply(mods, `[[`, "", "proteina")))
  cds_ss <- vapply(mods, function(m) substr(m$cds_sec, 1, 3 * nchar(m$proteina)), "")
  aln_cds <- alinear_codones(cds_ss, aln_prot)
  if (modo == "gen") {
    gen_sec <- vapply(mods, function(m) substr(m$sec, max(1, m$tss - promotor), m$tes), "")
    cat("  Alineando genes completos (puede tardar unos segundos)...\n")
    aln_gen <- alinear(DNAStringSet(gen_sec))
    ids <- list("Gen completo" = identidad(aln_gen), "CDS" = identidad(aln_cds),
                "Proteína" = identidad(aln_prot))
  } else {
    cat(sprintf("\n  El gen mide hasta %s nt: alinear los intrones completos sería muy lento.\n", format(largo_max, big.mark = ".")),
        " Las Figuras 1 y 2 muestran la REGIÓN CODIFICANTE (exones unidos).\n",
        " Para ver un exón con sus bordes intrónicos usa alinear_exon(datos, exon = n).\n\n")
    aln_gen <- NULL
    ids <- list("CDS" = identidad(aln_cds), "Proteína" = identidad(aln_prot))
  }
  datos$aln <- list(gen = aln_gen, cds = aln_cds, prot = aln_prot, promotor = promotor, modo = modo)
  datos$identidad <- ids

  # Matriz y anotación de la referencia sobre las columnas
  M <- matriz_aln(if (modo == "gen") aln_gen else aln_cds)
  pos_ref <- columnas_a_posicion(M[ref, ])
  hay <- which(!is.na(pos_ref))
  if (modo == "gen") {
    off <- max(1, mr$tss - promotor) - 1
    reg_col <- ifelse(is.na(pos_ref), NA, mr$region[off + pos_ref])
  } else {
    off <- 0
    exon_de_cds <- findInterval(mr$mapa_cds, mr$exones$inicio)
    ex <- exon_de_cds[pmin(pos_ref, length(exon_de_cds))]
    reg_col <- ifelse(is.na(pos_ref), NA, ifelse(ex %% 2 == 1, "exón impar", "exón par"))
  }
  for (k in seq_along(reg_col)) if (is.na(reg_col[k]))
    reg_col[k] <- if (k < min(hay) || k > max(hay)) "intergénico" else reg_col[k - 1]
  datos$reg_col <- reg_col
  datos$aln$M <- M; datos$aln$pos_ref <- pos_ref; datos$aln$off <- off

  nombres <- nombre_esp(esp)
  df <- data.frame(col = rep(seq_len(ncol(M)), times = nrow(M)),
                   especie = factor(rep(nombres, each = ncol(M)), levels = rev(nombres)),
                   base = as.vector(t(M)))
  refbase <- rep(M[ref, ], times = nrow(M))
  df$dif <- ifelse(df$base == "-", "-", ifelse(df$base == refbase, "igual", df$base))
  datos$aln$df <- df
  colores_reg <- c(COL_REG, "exón impar" = "#1F4E79", "exón par" = "#6FA3C8")
  anot <- data.frame(col = seq_len(ncol(M)), y = "Estructura",
                     reg = factor(reg_col, levels = names(colores_reg)))

  p_anot <- ggplot(anot, aes(col, y, fill = reg)) + geom_raster() +
    scale_fill_manual(values = colores_reg, drop = TRUE, name = paste0("Región (", nombre_esp(ref), ")")) +
    scale_x_continuous(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
    tema_base(11) + theme(axis.text.x = element_blank(), legend.position = "top")
  p_todo <- ggplot(df, aes(col, especie, fill = base)) + geom_raster() +
    scale_fill_manual(values = COL_NT, name = "Nucleótido") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, subtitle = "A. Todos los nucleótidos coloreados") +
    tema_base(11) + theme(axis.text.x = element_blank())
  p_dif <- ggplot(df, aes(col, especie, fill = dif)) + geom_raster() +
    scale_fill_manual(values = c(COL_NT, igual = "#F4F4F4"), name = NULL,
                      breaks = c("A", "C", "G", "T", "-"),
                      labels = c("A", "C", "G", "T", "hueco (indel)")) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Posición en el alineamiento (columna)", y = NULL,
         subtitle = paste0("B. Solo se colorean las diferencias con ", nombre_esp(ref), " (gris = igual)")) +
    tema_base(11)
  p1 <- (p_anot / p_todo / p_dif) + plot_layout(heights = c(0.6, length(esp), length(esp))) +
    plot_annotation(
      title = paste0("Figura 1. Alineamiento del gen ", datos$gen, " en ", length(esp), " especies"),
      subtitle = if (modo == "gen") "Promotor + exones + intrones. Los exones forman bandas conservadas; los intrones acumulan diferencias y huecos."
                 else "Región codificante (exones unidos; los tonos de azul alternan entre exones consecutivos).",
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  mostrar(p1)
  guardar_figura(p1, "fig1_mosaico_alineamiento", 15, 3 + 1.1 * length(esp))

  # Figura 1b: acercamiento automático a un límite exón-intrón (o exón-exón)
  if (modo == "gen") {
    destino <- if (nrow(mr$intrones)) mr$intrones$inicio[1] - off else mr$atg - off
  } else {
    destino <- if (nrow(mr$exones) > 1) sum(mr$mapa_cds <= mr$exones$fin[1]) + 1 else 1
  }
  centro <- which(pos_ref == destino)[1]
  pz <- dibujar_zoom(datos, centro, ancho = 30,
    titulo = if (modo == "gen") "Figura 1b. Acercamiento: final del exón 1 e inicio del intrón 1"
             else "Figura 1b. Acercamiento: unión entre el exón 1 y el exón 2 (región codificante)",
    subtitulo = "Observa qué se conserva justo en el límite (¿empieza el intrón por GT en todas las especies?)")
  mostrar(pz)
  guardar_figura(pz, "fig1b_acercamiento_limite", 15, 2.5 + 0.55 * length(esp))

  # Figura 2: perfil de identidad con la referencia
  otros <- setdiff(esp, ref)
  perf <- do.call(rbind, lapply(otros, function(e) {
    a <- M[ref, ]; b <- M[e, ]
    igual <- as.numeric(a == b & a != "-")
    igual[a == "-" & b == "-"] <- NA
    suav <- as.numeric(stats::filter(igual, rep(1 / ventana, ventana), sides = 2))
    data.frame(col = seq_along(a), identidad = 100 * suav, especie = nombre_esp(e))
  }))
  perf$especie <- factor(perf$especie, levels = nombre_esp(otros))
  es_exon <- reg_col %in% c("CDS", "5′UTR", "3′UTR", "exón impar")
  rr <- rle(es_exon); fin <- cumsum(rr$lengths); ini <- fin - rr$lengths + 1
  rect <- data.frame(xmin = ini[rr$values], xmax = fin[rr$values])
  p2 <- ggplot(perf, aes(col, identidad, colour = especie)) +
    geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#DCE9F2", alpha = 0.8) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = c("#7B3294", "#E66101", "#008837", "#2C7BB6", "#A6611A", "#5E3C99")[seq_along(otros)], name = NULL) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Posición en el alineamiento (columna)",
         y = paste0("% identidad con ", nombre_esp(ref)),
         title = paste0("Figura 2. Conservación a lo largo del gen ", datos$gen),
         subtitle = paste0("Ventana móvil de ", ventana, " columnas. Franjas azules = ",
                           if (modo == "gen") paste0("exones de ", nombre_esp(ref), ".") else "exones impares (1, 3, 5...).")) +
    tema_base(12)
  mostrar(p2)
  guardar_figura(p2, "fig2_perfil_conservacion", 13, 5)

  # Figura 3: matrices de identidad
  hm <- do.call(rbind, lapply(names(ids), function(t) {
    I <- ids[[t]][esp, esp]
    data.frame(f = rep(nombres, times = length(esp)), c = rep(nombres, each = length(esp)),
               v = as.vector(I), tipo = t)
  }))
  hm$f <- factor(hm$f, levels = rev(nombres)); hm$c <- factor(hm$c, levels = nombres)
  hm$tipo <- factor(hm$tipo, levels = names(ids))
  p3 <- ggplot(hm, aes(c, f, fill = v)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f", v), colour = v > 75), size = 4.2) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1B2A33"), guide = "none") +
    scale_fill_gradient(low = "#F7F1E3", high = "#1F4E79", limits = c(0, 100), name = "% identidad") +
    facet_wrap(~tipo) + labs(x = NULL, y = NULL,
      title = "Figura 3. Porcentaje de identidad entre especies",
      subtitle = "¿Qué se conserva más: el gen completo, la región codificante o la proteína?") +
    tema_base(12) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  mostrar(p3)
  guardar_figura(p3, "fig3_matriz_identidad", 4 + 3.4 * length(ids), 5.2)

  # Resumen numérico
  res <- data.frame(Especie = nombre_esp(otros),
                    `Divergencia (Ma)` = ESPECIES_INFO$divergencia_ma[match(otros, ESPECIES_INFO$clave)],
                    `Gen %` = if (modo == "gen") round(ids[["Gen completo"]][ref, otros], 1) else NA,
                    `Exones %` = NA, `Intrones %` = NA,
                    `CDS %` = round(ids[["CDS"]][ref, otros], 1),
                    `Proteína %` = round(ids[["Proteína"]][ref, otros], 1), check.names = FALSE)
  if (modo == "gen") {
    es_intron <- reg_col == "intrón"
    for (k in seq_along(otros)) {
      a <- M[ref, ]; b <- M[otros[k], ]
      f <- function(sel) { cols <- sel & (a != "-" | b != "-"); round(100 * mean(a[cols] == b[cols]), 1) }
      res$`Exones %`[k] <- f(es_exon)
      res$`Intrones %`[k] <- if (any(es_intron)) f(es_intron) else NA
    }
  }
  datos$resumen_identidad <- res
  cat("\n  Identidad con", nombre_esp(ref), "(%):\n")
  print(res, row.names = FALSE)
  write.csv(res, ruta("resultados", paste0(datos$gen, "_identidad.csv")), row.names = FALSE)
  cat("\n")
  preguntas(
    "¿Disminuye la identidad  a medida que aumenta el tiempo de divergencia? Ordena las especies.",
    "Compara la identidad de exones e intrones. ¿Por qué los intrones cambian más?",
    "¿Por qué la proteína puede ser más parecida que la secuencia de ADN que la codifica?",
    "En la Figura 1b, ¿qué nucleótidos se conservan en el límite exón|intrón?",
    "En la Figura 2, ¿qué exón es el más conservado? ¿Qué parte de la proteína codifica? (ver Figura 4)",
    "Si el reloj molecular fuera perfecto, ¿Caería la identidad en proporción al tiempo? Grafica % identidad CDS vs. millones de años con los datos de la tabla.")
  preguntas_busqueda(
    paste0("Usa acercamiento(datos) y escribe la posición de la zona MENOS conservada de la Figura 2. ",
           "¿Está en un intrón, en un UTR o en la CDS? Describe qué ves."),
    paste0("En NCBI Gene (", url_ncbi_gene(datos$gen), ") abre la pestaña 'Orthologs': ",
           "¿cuántas especies tienen este gen? ¿Aparecen especies en las que no se encontró?"),
    "Busca el término 'selección purificadora' (purifying selection): ¿cómo explica la diferencia entre exones e intrones?")
  invisible(datos)
}

# Dibuja un acercamiento con letras alrededor de una columna del alineamiento.
dibujar_zoom <- function(datos, centro, ancho = 30, marca = NULL, titulo = "", subtitulo = "") {
  M <- datos$aln$M; df <- datos$aln$df; reg_col <- datos$reg_col
  esp <- orden_especies(datos)
  colores_reg <- c(COL_REG, "exón impar" = "#1F4E79", "exón par" = "#6FA3C8")
  c1 <- max(1, centro - ancho); c2 <- min(ncol(M), centro + ancho)
  dz <- df[df$col >= c1 & df$col <= c2, ]
  az <- data.frame(col = c1:c2, reg = factor(reg_col[c1:c2], levels = names(colores_reg)))
  p <- ggplot() +
    geom_tile(data = dz, aes(col, especie, fill = base), colour = "white") +
    geom_text(data = dz, aes(col, especie, label = base), size = 3.6, colour = "white", fontface = "bold") +
    geom_tile(data = az, aes(col, length(esp) + 1, fill = reg), height = 0.5) +
    scale_fill_manual(values = c(COL_NT, colores_reg),
                      breaks = intersect(c("promotor", "5′UTR", "CDS", "intrón", "3′UTR", "exón impar", "exón par"),
                                         as.character(az$reg)), name = "Región") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
    labs(x = "Columna del alineamiento", y = NULL, title = titulo, subtitle = subtitulo) +
    tema_base(12)
  if (!is.null(marca))
    p <- p + annotate("rect", xmin = marca - 0.5, xmax = marca + 0.5, ymin = 0.4,
                      ymax = length(esp) + 0.6, fill = NA, colour = "#1B2A33", linewidth = 1.2)
  p
}

# ---------- FIGURA 1c: acercamiento a una posición o a un rango ---------------
# posicion: un número relativo al TSS (350, -30), una posición en la CDS ("c.20"),
# o un RANGO de dos posiciones: c(300, 400), c("c.20", "c.80"), "300:400" o "c.20-c.80".

# Convierte una posición escrita por el estudiante en coordenadas del gen.
interpretar_posicion <- function(datos, x) {
  mr <- datos$modelos[[datos$ref]]
  x <- gsub("\\s", "", as.character(x))
  es_c <- grepl("^c\\.", x)
  n <- suppressWarnings(as.integer(sub("^c\\.", "", x)))
  if (is.na(n)) stop("No entendí la posición '", x, "'. Usa un número (350, -30) o c.N (c.20).")
  if (es_c) {
    if (n < 1 || n > length(mr$mapa_cds)) stop("La CDS va de c.1 a c.", length(mr$mapa_cds))
    g <- mr$mapa_cds[n]; idx <- n
  } else { g <- mr$tss + n; idx <- match(g, mr$mapa_cds) }
  if (datos$aln$modo == "gen") {
    if (g < mr$tss - datos$aln$promotor || g > mr$tes)
      stop(sprintf("La posición %s está fuera de la región alineada (de %d a +%d respecto al TSS).",
                   x, -datos$aln$promotor, mr$tes - mr$tss))
    objetivo <- g - datos$aln$off
  } else {
    if (is.na(idx)) stop("En genes largos solo se alinea la región codificante: usa posiciones c.N (entre c.1 y c.", length(mr$mapa_cds), ").")
    objetivo <- idx
  }
  col <- which(datos$aln$pos_ref == objetivo)[1]
  region <- mr$region[g]
  donde <- if (region == "intrón") paste("intrón", which(mr$intrones$inicio <= g & mr$intrones$fin >= g))
           else if (region %in% c("CDS", "5′UTR", "3′UTR")) paste0(region, ", exón ", which(mr$exones$inicio <= g & mr$exones$fin >= g))
           else region
  list(texto = x, g = g, idx = idx, col = col, donde = donde,
       etiqueta = sprintf("%+d%s (%s)", g - mr$tss, if (!is.na(idx)) paste0(" = c.", idx) else "", donde))
}

acercamiento <- function(datos, posicion = NULL, ancho = 30) {
  if (is.null(datos$aln)) stop("Primero ejecuta el MÓDULO 1 (alineamiento).")
  mr <- datos$modelos[[datos$ref]]
  if (is.null(posicion) || identical(posicion, "")) {
    if (!interactive()) stop("Indica la posición, por ejemplo acercamiento(datos, 350), acercamiento(datos, 'c.20') o acercamiento(datos, c(300, 400)).")
    cat("\n  Escribe la posición o el rango que quieres ver en el gen de", nombre_esp(datos$ref), "\n",
        "  - Una posición relativa al TSS (0 = inicio del ARNm; negativos = promotor). Ej.: 350  o  -30\n",
        "  - Una posición en la región codificante con 'c.'. Ej.: c.20\n",
        "  - Un RANGO con dos posiciones separadas por ':'. Ej.: 300:400  o  c.20:c.80\n",
        sprintf("  (El transcrito va de 0 a +%d; la CDS de c.1 a c.%d)\n", mr$tes - mr$tss, nchar(mr$cds_sec)))
    posicion <- readline("  Posición o rango: ")
  }
  # ¿Una posición o un rango?
  if (length(posicion) == 1) {
    txt <- gsub("\\s", "", as.character(posicion))
    partes <- strsplit(txt, ":|(?<=[0-9])-(?=c\\.|[0-9])", perl = TRUE)[[1]]
    posicion <- partes[nzchar(partes)]
  }
  if (length(posicion) > 2) stop("Escribe una posición o un rango de DOS posiciones, p. ej. c(300, 400).")

  if (length(posicion) == 1) {
    p <- interpretar_posicion(datos, posicion)
    centro <- p$col; c1 <- max(1, centro - ancho); c2 <- min(ncol(datos$aln$M), centro + ancho)
    marcas <- centro
    info <- paste0("Posición ", p$etiqueta)
    nombre <- paste0("fig1c_acercamiento_", gsub("[^A-Za-z0-9-]", "", p$texto))
    idx_cod <- p$idx
  } else {
    p1 <- interpretar_posicion(datos, posicion[1]); p2 <- interpretar_posicion(datos, posicion[2])
    if (p1$col > p2$col) { tmp <- p1; p1 <- p2; p2 <- tmp }
    c1 <- p1$col; c2 <- p2$col
    if (c2 - c1 > 600) { c2 <- c1 + 600; message("El rango es muy largo: se muestran las primeras 600 columnas.") }
    marcas <- c(c1, c2)
    info <- paste0("Rango de ", p1$etiqueta, " a ", p2$etiqueta, " · ", c2 - c1 + 1, " columnas del alineamiento")
    nombre <- paste0("fig1c_rango_", gsub("[^A-Za-z0-9-]", "", p1$texto), "_a_", gsub("[^A-Za-z0-9-]", "", p2$texto))
    idx_cod <- if (!is.na(p1$idx) && !is.na(p2$idx)) c(p1$idx, p2$idx) else NA
  }

  # Aminoácidos de cada especie en el codón (posición única) o resumen del tramo (rango)
  esp <- orden_especies(datos)
  P <- matriz_aln(datos$aln$prot[esp]); pp <- columnas_a_posicion(P[datos$ref, ])
  if (length(idx_cod) == 1 && !is.na(idx_cod)) {
    cod_n <- (idx_cod - 1) %/% 3 + 1
    col_p <- which(pp == cod_n)[1]
    info <- paste0(info, "\nCodón ", cod_n, " -> ", paste0(nombre_esp(esp), ": ", aa3(P[, col_p]), collapse = " · "))
  } else if (length(idx_cod) == 2) {
    k1 <- (idx_cod[1] - 1) %/% 3 + 1; k2 <- (idx_cod[2] - 1) %/% 3 + 1
    cols <- which(pp >= k1 & pp <= k2)
    difs <- vapply(esp[-1], function(e) sum(P[e, cols] != P[datos$ref, cols]), 1)
    info <- paste0(info, "\nCodones ", k1, "-", k2, " · aminoácidos distintos a ", nombre_esp(datos$ref), ": ",
                   paste0(nombre_esp(esp[-1]), " ", difs, collapse = " · "))
  }
  cat("\n ", gsub("\n", "\n  ", info), "\n\n")

  if (length(marcas) == 1) {
    p <- dibujar_zoom(datos, marcas, ancho, marca = marcas,
                      titulo = paste0("Figura 1c. Acercamiento a la posición elegida en ", datos$gen),
                      subtitulo = paste0(info, "\nEl recuadro negro marca la posición elegida."))
    alto <- 3 + 0.55 * length(esp)
  } else {
    p <- dibujar_rango(datos, c1, c2,
                       titulo = paste0("Figura 1c. Tramo elegido del gen ", datos$gen),
                       subtitulo = paste0(info, "\nLos recuadros negros marcan el inicio y el final del rango."))
    nfil <- ceiling((c2 - c1 + 1) / 100)
    alto <- 2.2 + nfil * (0.42 * length(esp) + 0.9)
  }
  mostrar(p)
  guardar_figura(p, nombre, 15, alto)
  preguntas("¿En qué región del gen está la posición o el tramo? ¿Es igual en todas las especies?",
            "¿Las diferencias que ves son sustituciones o huecos (indels)? ¿Se concentran en alguna parte del tramo?")
  invisible(datos)
}

# Dibuja un tramo del alineamiento (columnas c1 a c2) en filas de 100 columnas.
dibujar_rango <- function(datos, c1, c2, titulo = "", subtitulo = "", por_fila = 100) {
  df <- datos$aln$df; reg_col <- datos$reg_col
  esp <- orden_especies(datos); ne <- length(esp)
  colores_reg <- c(COL_REG, "exón impar" = "#1F4E79", "exón par" = "#6FA3C8")
  d <- df[df$col >= c1 & df$col <= c2, ]
  d$x <- (d$col - c1) %% por_fila + 1
  d$fila <- (d$col - c1) %/% por_fila
  a <- data.frame(col = c1:c2, reg = factor(reg_col[c1:c2], levels = names(colores_reg)))
  a$x <- (a$col - c1) %% por_fila + 1; a$fila <- (a$col - c1) %/% por_fila
  etiqueta_fila <- function(f) paste0("Columnas ", c1 + f * por_fila, "-", pmin(c2, c1 + (f + 1) * por_fila - 1))
  d$tramo <- factor(etiqueta_fila(d$fila), levels = etiqueta_fila(sort(unique(d$fila))))
  a$tramo <- factor(etiqueta_fila(a$fila), levels = levels(d$tramo))
  mk <- data.frame(col = c(c1, c2)); mk$x <- (mk$col - c1) %% por_fila + 1
  mk$tramo <- factor(etiqueta_fila((mk$col - c1) %/% por_fila), levels = levels(d$tramo))
  ggplot() +
    geom_tile(data = d, aes(x, especie, fill = base), colour = "white", linewidth = 0.2) +
    geom_text(data = d, aes(x, especie, label = base), size = 2.4, colour = "white", fontface = "bold") +
    geom_tile(data = a, aes(x, ne + 0.95, fill = reg), height = 0.45) +
    geom_tile(data = mk, aes(x, (ne + 1) / 2), width = 1, height = ne, fill = NA, colour = "#1B2A33", linewidth = 0.9) +
    scale_fill_manual(values = c(COL_NT, colores_reg),
                      breaks = intersect(c("promotor", "5′UTR", "CDS", "intrón", "3′UTR", "exón impar", "exón par"),
                                         as.character(a$reg)), name = "Región") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, por_fila + 0.5)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1))) +
    facet_wrap(~tramo, ncol = 1) +
    labs(x = NULL, y = NULL, title = titulo, subtitle = subtitulo) +
    tema_base(11) + theme(axis.text.x = element_blank(), strip.text = element_text(hjust = 0))
}

# ---------- MÓDULO 2: estructura del gen -------------------------------------
modulo2_estructura <- function(datos, promotor = 200) {
  caja("MÓDULO 2 · Estructura del gen: promotor, exones, intrones y señales",
       "TSS = sitio de inicio de la transcripción (inicio del ARNm).",
       "ATG = codón de inicio de la traducción. ¡No es lo mismo que el TSS!",
       "Entre el TSS y el ATG está el 5′UTR; después del codón de parada, el 3′UTR.",
       "Cada intrón empieza en GT (sitio donador) y termina en AG (sitio aceptor).")
  esp <- orden_especies(datos)
  mods <- datos$modelos[esp]
  nombres <- nombre_esp(esp)

  # Tabla de estructura
  tab <- data.frame(
    Especie = nombres,
    `Gen TSS-TES (nt)` = vapply(mods, function(m) m$tes - m$tss + 1, 1),
    Exones = vapply(mods, function(m) nrow(m$exones), 1),
    `Exones (nt)` = vapply(mods, function(m) length(m$mapa_arn), 1),
    `Intrones (nt)` = vapply(mods, function(m) sum(m$intrones$fin - m$intrones$inicio + 1), 1),
    `5'UTR (nt)` = vapply(mods, function(m) sum(m$mapa_arn < m$atg), 1),
    `CDS (nt)` = vapply(mods, function(m) length(m$mapa_cds), 1),
    `3'UTR (nt)` = vapply(mods, function(m) sum(m$mapa_arn > m$fin_cds), 1),
    check.names = FALSE)
  tab$`% intrónico` <- round(100 * tab$`Intrones (nt)` / tab$`Gen TSS-TES (nt)`, 1)
  datos$tabla_estructura <- tab
  print(tab, row.names = FALSE)
  write.csv(tab, ruta("resultados", paste0(datos$gen, "_estructura.csv")), row.names = FALSE)

  # Figura 4: mapas de los genes, en coordenadas relativas al TSS
  bloques <- do.call(rbind, lapply(seq_along(esp), function(i) {
    m <- mods[[i]]
    do.call(rbind, lapply(seq_len(nrow(m$exones)), function(k) {
      e <- m$exones[k, ]
      p <- e$inicio:e$fin
      codif <- p %in% m$mapa_cds
      rr <- rle(codif); f <- cumsum(rr$lengths); s <- f - rr$lengths + 1
      data.frame(y = i * 1.7, xmin = p[s] - m$tss - 0.5, xmax = p[f] - m$tss + 0.5,
                 cds = rr$values, exon = k)
    }))
  }))
  bloques$alto <- ifelse(bloques$cds, 0.32, 0.16)
  ESP <- 1.7   # separación vertical entre especies
  lineas <- data.frame(y = seq_along(esp) * ESP,
                       x0 = -promotor, x1 = vapply(mods, function(m) m$tes - m$tss, 1))
  marcas <- do.call(rbind, lapply(seq_along(esp), function(i) {
    m <- mods[[i]]
    d <- data.frame(y = i * ESP, x = c(m$atg - m$tss, m$fin_cds - 2 - m$tss),
                    lab = c("ATG", "Stop"))
    if (!is.na(m$poliA)) d <- rbind(d, data.frame(y = i * ESP, x = m$poliA - m$tss, lab = "AATAAA"))
    d
  }))
  etiq_ex <- aggregate(cbind(xmin, xmax) ~ y + exon, data = bloques, FUN = function(v) v[1])
  etiq_ex$xmin <- aggregate(xmin ~ y + exon, data = bloques, FUN = min)$xmin
  etiq_ex$xmax <- aggregate(xmax ~ y + exon, data = bloques, FUN = max)$xmax
  intr <- do.call(rbind, lapply(seq_along(esp), function(i) {
    m <- mods[[i]]
    if (!nrow(m$intrones)) return(NULL)
    data.frame(y = i * ESP, x = (m$intrones$inicio + m$intrones$fin) / 2 - m$tss,
               lab = paste0(m$intrones$fin - m$intrones$inicio + 1, " nt"))
  }))
  escala <- max(lineas$x1)
  p4 <- ggplot() +
    geom_segment(data = lineas, aes(x = x0, xend = x1, y = y, yend = y), colour = "#7A7A7A", linewidth = 0.6) +
    geom_segment(data = lineas, aes(x = -promotor, xend = 0, y = y, yend = y), colour = "#B5A98A", linewidth = 2) +
    geom_rect(data = bloques, aes(xmin = xmin, xmax = xmax, ymin = y - alto, ymax = y + alto,
                                  fill = ifelse(cds, "Región codificante (CDS)", "UTR (exón no traducido)"))) +
    geom_text(data = etiq_ex, aes(x = (xmin + xmax) / 2, y = y + 0.45, label = paste0("E", exon)), size = 3.4, fontface = "bold") +
    geom_segment(data = data.frame(y = seq_along(esp) * ESP), aes(x = 0, xend = 0, y = y + 0.05, yend = y + 0.62), linewidth = 0.6) +
    geom_segment(data = data.frame(y = seq_along(esp) * ESP), aes(x = 0, xend = escala * 0.04, y = y + 0.62, yend = y + 0.62),
                 linewidth = 0.6, arrow = arrow(length = unit(0.12, "cm"))) +
    geom_text(data = data.frame(y = seq_along(esp) * ESP), aes(x = 0, y = y + 0.78, label = "TSS"), size = 3.2, hjust = 0.2) +
    geom_point(data = marcas, aes(x = x, y = y - 0.5, shape = lab, colour = lab), size = 3) +
    scale_shape_manual(values = c(ATG = 17, Stop = 15, AATAAA = 18), name = "Señal") +
    scale_colour_manual(values = c(ATG = "#2E9E44", Stop = "#D62728", AATAAA = "#7B3294"), name = "Señal") +
    scale_fill_manual(values = c("Región codificante (CDS)" = "#1F4E79", "UTR (exón no traducido)" = "#8DB8D0"), name = NULL) +
    scale_y_continuous(breaks = seq_along(esp) * ESP, labels = nombres, limits = c(ESP - 0.7, length(esp) * ESP + 0.95)) +
    labs(x = "Posición relativa al TSS (nt)", y = NULL,
         title = paste0("Figura 4. Estructura del gen ", datos$gen, " en cada especie"),
         subtitle = paste0("Barra beige = promotor (", promotor, " nt antes del TSS). Línea gris = intrones (con su longitud). Cajas = exones.")) +
    tema_base(12) + theme(panel.grid.major.x = element_line(colour = "#EEEEEE"))
  if (!is.null(intr))
    p4 <- p4 + geom_text(data = intr, aes(x = x, y = y - 0.22, label = lab), size = 3, colour = "#555555")
  mostrar(p4)
  guardar_figura(p4, "fig4_estructura_gen", 14, 1.8 + 1.5 * length(esp))

  # Figura 5: límites exón-intrón (sitios donador y aceptor) en todas las especies
  nint <- min(2, max(vapply(mods, function(m) nrow(m$intrones), 1)))
  if (nint == 0) {
    cat("\n  Este gen no tiene intrones en su región codificante: no hay límites exón-intrón.\n",
        " Pregunta: ¿qué ventaja o desventaja tiene un gen sin intrones?\n\n")
  } else {
    filas <- list()
    for (i in seq_along(esp)) {
      m <- mods[[i]]
      for (k in seq_len(min(nint, nrow(m$intrones)))) {
        a <- m$intrones$inicio[k]; b <- m$intrones$fin[k]
        don <- (a - 8):(a + 9); acc <- (b - 9):(b + 8)
        for (tipo in c("donador", "aceptor")) {
          pos <- if (tipo == "donador") don else acc
          pos <- pos[pos >= 1 & pos <= m$L]
          esintr <- pos >= a & pos <= b
          clave <- if (tipo == "donador") pos %in% c(a, a + 1) else pos %in% c(b - 1, b)
          filas[[length(filas) + 1]] <- data.frame(
            especie = nombres[i], panel = paste0("Intrón ", k, " · sitio ", tipo,
                                                 if (tipo == "donador") " (exón | intrón)" else " (intrón | exón)"),
            x = seq_along(pos), letra = ifelse(esintr, tolower(m$nt[pos]), m$nt[pos]),
            zona = ifelse(clave, "GT / AG", ifelse(esintr, "Intrón", "Exón")))
        }
      }
    }
    dj <- do.call(rbind, filas)
    dj$especie <- factor(dj$especie, levels = rev(nombres))
    dj$panel <- factor(dj$panel, levels = unique(dj$panel))
    p5 <- ggplot(dj, aes(x, especie)) +
      geom_tile(aes(fill = zona), colour = "white", linewidth = 0.6) +
      geom_text(aes(label = letra, fontface = ifelse(zona == "GT / AG", "bold", "plain")), size = 4.2) +
      scale_fill_manual(values = c("Exón" = "#BFD9E8", "Intrón" = "#F4F1EA", "GT / AG" = "#F6C343"), name = NULL) +
      facet_wrap(~panel, ncol = 2, scales = "free_x") +
      labs(x = NULL, y = NULL, title = "Figura 5. Límites exón-intrón en todas las especies",
           subtitle = "Exón en MAYÚSCULAS (azul), intrón en minúsculas. Resaltado: los dinucleótidos de los sitios de splicing.") +
      tema_base(12) + theme(axis.text.x = element_blank(), strip.text = element_text(face = "bold"))
    mostrar(p5)
    guardar_figura(p5, "fig5_limites_exon_intron", 14, 1.8 + 0.5 * length(esp) * ceiling(nint))
  }

  # Figura 5b: secuencia anotada de la especie de referencia
  mr <- mods[[datos$ref]]
  a <- max(1, mr$tss - 60); b <- min(mr$L, a + 899)
  pos <- a:b
  ds <- data.frame(pos = pos, letra = mr$nt[pos], reg = mr$region[pos])
  ds$letra <- ifelse(ds$reg == "intrón" | ds$reg == "promotor" | ds$reg == "intergénico", tolower(ds$letra), ds$letra)
  ds$fila <- (seq_along(pos) - 1) %/% 60
  ds$col <- (seq_along(pos) - 1) %% 60 + 1
  ds$reg <- factor(ds$reg, levels = names(COL_REG))
  marcas_s <- c(TSS = mr$tss, ATG = mr$atg)
  ds$marca <- ifelse(ds$pos %in% (mr$atg + 0:2), "ATG", ifelse(ds$pos == mr$tss, "TSS", ""))
  etiq_fila <- data.frame(fila = unique(ds$fila), pos = tapply(ds$pos, ds$fila, min) - mr$tss)
  p5b <- ggplot(ds, aes(col, -fila)) +
    geom_tile(aes(fill = reg), colour = "white", linewidth = 0.3) +
    geom_tile(data = ds[ds$marca != "", ], fill = NA, colour = "#D62728", linewidth = 0.9) +
    geom_text(aes(label = letra, colour = reg == "CDS"), size = 3.1, family = "mono") +
    geom_text(data = etiq_fila, aes(x = -2, y = -fila, label = sprintf("%+d", pos)), size = 3, hjust = 1, colour = "#555555") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1B2A33"), guide = "none") +
    scale_fill_manual(values = COL_REG, drop = TRUE, name = "Región") +
    scale_x_continuous(limits = c(-6, 61), expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 5b. Secuencia anotada de ", datos$gen, " (", nombre_esp(datos$ref), "), desde 60 nt antes del TSS"),
         subtitle = "Números = posición relativa al TSS. Recuadro rojo = TSS y codón ATG. Exones en MAYÚSCULAS, intrones y promotor en minúsculas.") +
    tema_base(11) + theme(axis.text = element_blank())
  mostrar(p5b)
  guardar_figura(p5b, "fig5b_secuencia_anotada", 14, 1.5 + 0.3 * (max(ds$fila) + 1))

  # Archivos de texto con la secuencia anotada de cada especie
  for (e in esp) {
    m <- mods[[e]]
    s <- m$nt
    minus <- m$region %in% c("intrón", "promotor", "intergénico")
    s[minus] <- tolower(s[minus])
    sub_s <- s[max(1, m$tss - promotor):m$tes]
    lineas_txt <- vapply(split(sub_s, ceiling(seq_along(sub_s) / 60)), paste, "", collapse = "")
    writeLines(c(paste0("> ", datos$gen, " ", nombre_esp(e), " | ", m$accesion,
                        " | promotor (", promotor, " nt) + gen. MAYÚSCULAS = exones; minúsculas = promotor e intrones"),
                 lineas_txt),
               ruta("resultados", paste0(datos$gen, "_", e, "_secuencia_anotada.txt")))
  }
  cat("\n  Secuencias anotadas guardadas en resultados/*_secuencia_anotada.txt\n\n")
  preguntas(
    "¿Cuántos exones tiene el gen en cada especie? ¿Cambia el número o solo la longitud de los intrones?",
    "¿Qué porcentaje del gen corresponde a intrones? ¿Qué especie tiene el gen más compacto?",
    "¿A cuántos nucleótidos del TSS está el ATG? ¿Cómo se llama esa región?",
    "¿Se cumple la regla GT-AG en todos los intrones y especies (Figura 5)?",
    "¿Coincide la posición de los intrones dentro de la proteína entre especies? ¿Qué sugiere sobre su antigüedad?",
    "Compara la longitud de los UTR entre especies: ¿qué región del ARNm varía más, los UTR o la CDS?")
  preguntas_busqueda(
    "Ejecuta buscar_secuencia(datos, 'GTAAGT') y buscar_secuencia(datos, 'GTRAGT'): ¿cuántos intrones tienen el consenso completo del sitio donador?",
    "Ejecuta perfil_gc(datos): ¿el promotor o el exón 1 cumplen el criterio de isla CpG? ¿Por qué importa para la metilación?",
    paste0("En NCBI Datasets (", url_genoma, ") busca el tamaño del genoma de cada especie. ",
           "¿Las especies con genoma más grande tienen intrones más largos en este gen?"),
    paste0("En NCBI Gene (", url_ncbi_gene(datos$gen), ") revisa cuántas variantes de transcrito (isoformas) tiene el gen humano."))
  invisible(datos)
}

# ---------- Utilidades para dibujar secuencias con letras --------------------

# Divide una secuencia en filas de 'ancho' posiciones, sin partir codones
# (si se da el índice CDS de cada posición).
filas_por_codon <- function(idx_cds, ancho = 60) {
  n <- length(idx_cds)
  fila <- integer(n); col <- integer(n)
  f <- 0; c <- 0
  for (i in seq_len(n)) {
    inicio_codon <- is.na(idx_cds[i]) || (idx_cds[i] - 1) %% 3 == 0
    if (c >= ancho - 2 && inicio_codon && (n - i) >= 3) { f <- f + 1; c <- 0 }
    c <- c + 1
    fila[i] <- f; col[i] <- c
  }
  data.frame(fila = fila, col = col)
}

exon_info <- function(m, k) {
  k <- max(1, min(k, nrow(m$exones)))
  e <- m$exones[k, ]
  pos <- e$inicio:e$fin
  idx <- match(pos, m$mapa_cds)        # posición en la CDS (NA = UTR)
  list(k = k, pos = pos, idx = idx, sec = paste(m$nt[pos], collapse = ""))
}

# ---------- MÓDULO 3: transcripción ------------------------------------------
modulo3_transcripcion <- function(datos, exon = 1, especie = NULL, max_nt = 180) {
  especie <- especie %||% datos$ref
  m <- datos$modelos[[especie]]
  ex <- exon_info(m, exon)
  if (ex$k != exon) message("El gen solo tiene ", nrow(m$exones), " exones; se usa el exón ", ex$k)
  caja(paste0("MÓDULO 3 · Transcripción del exón ", ex$k, " de ", datos$gen, " (", nombre_esp(especie), ")"),
       "La ARN polimerasa lee la hebra MOLDE en dirección 3′->5′ y sintetiza el ARNm 5′->3′.",
       "El ARNm es complementario a la hebra molde e IGUAL a la hebra codificante, con U en lugar de T.",
       "Nota: en la célula se transcribe el gen completo (pre-ARNm con intrones);",
       "aquí mostramos solo el tramo del exón para ver el proceso.")
  n <- min(nchar(ex$sec), max_nt)
  cod <- letras(substr(ex$sec, 1, n))
  molde <- unname(c(A = "T", T = "A", G = "C", C = "G", N = "N")[cod])
  arn <- ifelse(cod == "T", "U", cod)
  cat("  Hebra codificante 5′-", paste(cod, collapse = ""), "-3′\n", sep = "")
  cat("  Hebra molde       3′-", paste(molde, collapse = ""), "-5′\n", sep = "")
  cat("  ARNm              5′-", paste(arn, collapse = ""), "-3′\n\n", sep = "")
  if (nchar(ex$sec) > max_nt) cat("  (Se muestran los primeros", max_nt, "de", nchar(ex$sec), "nt del exón)\n\n")

  pos <- seq_len(n)
  fila <- (pos - 1) %/% 60; col <- (pos - 1) %% 60 + 1
  d <- rbind(data.frame(col, fila, y = 3, l = cod, t = "codificante"),
             data.frame(col, fila, y = 2, l = molde, t = "molde"),
             data.frame(col, fila, y = 0.6, l = arn, t = "ARNm"))
  d$yy <- d$y - fila * 4.6
  et <- expand.grid(fila = unique(fila), t = c("Hebra codificante  5′→3′", "Hebra molde  3′→5′", "ARNm  5′→3′"))
  et$y <- c(3, 2, 0.6)[as.integer(et$t)] - et$fila * 4.6
  enlaces <- data.frame(x = col, y1 = 2.72 - fila * 4.6, y2 = 2.28 - fila * 4.6)
  p6 <- ggplot(d) +
    geom_segment(data = enlaces, aes(x = x, xend = x, y = y1, yend = y2), colour = "#BBBBBB") +
    geom_tile(aes(col, yy, fill = l), width = 0.92, height = 0.5) +
    geom_text(aes(col, yy, label = l), colour = "white", fontface = "bold", size = 3.4) +
    geom_text(data = et, aes(x = -0.5, y = y, label = t), hjust = 1, size = 3.4) +
    geom_segment(data = data.frame(f = unique(fila)), aes(x = 1, xend = 12, y = 1.3 - f * 4.6, yend = 1.3 - f * 4.6),
                 arrow = arrow(length = unit(0.15, "cm")), colour = "#7B3294") +
    geom_text(data = data.frame(f = unique(fila)), aes(x = 13, y = 1.3 - f * 4.6, label = "ARN polimerasa (síntesis 5′→3′)"),
              hjust = 0, size = 3, colour = "#7B3294") +
    scale_fill_manual(values = COL_NT, guide = "none") +
    scale_x_continuous(limits = c(-17, 61), expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 6. Transcripción del exón ", ex$k, " de ", datos$gen, " (", nombre_esp(especie), ")"),
         subtitle = "El ARNm se copia de la hebra molde: tiene la misma secuencia que la hebra codificante, pero con U (rojo) en lugar de T.") +
    tema_base(12) + theme(axis.text = element_blank())
  mostrar(p6)
  guardar_figura(p6, "fig6_transcripcion_exon", 14, 1.4 + 1.6 * (max(fila) + 1))
  datos$exon_ejemplo <- ex$k
  preguntas(
    "¿Por qué el ARNm tiene la misma secuencia que la hebra codificante y no que la molde?",
    "¿Qué enzima realiza la transcripción y en qué dirección sintetiza?",
    "Si en la célula se transcribe el gen completo, ¿qué proceso elimina los intrones del pre-ARNm?",
    "¿Qué modificaciones recibe el pre-ARNm en sus extremos 5′ y 3′ antes de salir del núcleo?")
  preguntas_busqueda(
    "Usa ver_region(datos, -40, 10): ¿qué secuencia encuentras unos 30 nt antes del TSS? ¿Qué proteína la reconoce?",
    "Ejecuta buscar_secuencia(datos, 'AATAAA'): ¿cuántas copias hay y cuál está a 10-30 nt del final del ARNm?")
  invisible(datos)
}

# ---------- MÓDULO 4: traducción (relación lineal ADN - ARNm - proteína) -----
modulo4_traduccion <- function(datos, exon = NULL, especie = NULL, max_nt = 180) {
  especie <- especie %||% datos$ref
  exon <- exon %||% datos$exon_ejemplo %||% 1
  m <- datos$modelos[[especie]]
  ex <- exon_info(m, exon)
  caja(paste0("MÓDULO 4 · Traducción del exón ", ex$k, " de ", datos$gen, " (", nombre_esp(especie), ")"),
       "El ribosoma lee el ARNm de 3 en 3 nucleótidos (codones), desde el AUG.",
       "Cada codón corresponde a un aminoácido: el orden de los codones en el gen",
       "es el mismo orden de los aminoácidos en la proteína (COLINEALIDAD).",
       "Un codón puede quedar PARTIDO entre dos exones: por eso hay que conocer",
       "el marco de lectura (fase) del exón antes de traducirlo.")
  n <- min(length(ex$pos), max_nt)
  pos <- ex$pos[seq_len(n)]; idx <- ex$idx[seq_len(n)]
  nt <- m$nt[pos]
  arn <- ifelse(nt == "T", "U", nt)
  lay <- filas_por_codon(idx)
  codon <- ifelse(is.na(idx), NA, (idx - 1) %/% 3 + 1)
  # Codones completos o partidos dentro del exón mostrado
  cds_full <- m$cds_sec
  n_cod <- floor(nchar(cds_full) / 3)
  cods <- sort(unique(na.omit(codon)))
  bloques <- do.call(rbind, lapply(cods, function(cn) {
    sel <- which(codon == cn)
    completo <- length(sel) == 3
    trip <- substr(cds_full, 3 * cn - 2, 3 * cn)
    aa <- if (cn <= n_cod) traducir(trip) else "?"
    data.frame(codon = cn, fila = lay$fila[sel[1]], x1 = min(lay$col[sel]), x2 = max(lay$col[sel]),
               completo = completo, triplete = trip, aa = aa)
  }))
  d <- data.frame(col = lay$col, fila = lay$fila, nt = nt, arn = arn, utr = is.na(idx))
  H <- 5.2
  d$y_adn <- 3 - d$fila * H; d$y_arn <- 2 - d$fila * H
  if (!is.null(bloques) && nrow(bloques)) {
    bloques$y <- 2 - bloques$fila * H
    bloques$yaa <- 0.5 - bloques$fila * H
    bloques$xm <- (bloques$x1 + bloques$x2) / 2
    bloques$lab <- ifelse(bloques$completo, aa3(bloques$aa), paste0("(", aa3(bloques$aa), ")"))
    bloques$par <- factor(bloques$codon %% 2)
  }
  filas_u <- unique(d$fila)
  et <- data.frame(fila = rep(filas_u, each = 3),
                   lab = rep(c("ADN (codificante)", "ARNm", "Proteína"), times = length(filas_u)),
                   y = rep(c(3, 2, 0.5), times = length(filas_u)) - rep(filas_u, each = 3) * H)
  p7 <- ggplot() +
    geom_text(data = et, aes(x = -0.5, y = y, label = lab), hjust = 1, size = 3.4, fontface = "bold")
  if (!is.null(bloques) && nrow(bloques)) {
    p7 <- p7 +
      geom_rect(data = bloques, aes(xmin = x1 - 0.48, xmax = x2 + 0.48, ymin = y - 0.36, ymax = y + 0.36,
                                    linetype = completo, fill = par), colour = "#1B2A33", alpha = 0.25) +
      geom_segment(data = bloques, aes(x = xm, xend = xm, y = y - 0.4, yend = yaa + 0.32), colour = "#888888") +
      geom_label(data = bloques, aes(x = xm, y = yaa, label = lab, alpha = completo),
                 fill = "#F6C343", size = 3, label.size = 0.2, label.padding = unit(0.12, "lines")) +
      geom_text(data = bloques[bloques$completo, ], aes(x = xm, y = yaa - 0.55, label = aa), size = 3, colour = "#555555")
  }
  p7 <- p7 +
    geom_tile(data = d, aes(col, y_adn, fill = ifelse(utr, "UTR", nt)), width = 0.92, height = 0.5) +
    geom_text(data = d, aes(col, y_adn, label = nt), colour = "white", size = 3.2, fontface = "bold") +
    geom_tile(data = d, aes(col, y_arn, fill = ifelse(utr, "UTR", arn)), width = 0.92, height = 0.5) +
    geom_text(data = d, aes(col, y_arn, label = arn), colour = "white", size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c(COL_NT, UTR = "#9E9E9E", `0` = "#FFFFFF", `1` = "#BFD9E8"), guide = "none") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
    scale_x_continuous(limits = c(-12, 62), expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 7. Del ADN a la proteína: exón ", ex$k, " de ", datos$gen, " (", nombre_esp(especie), ")"),
         subtitle = "Gris = UTR (no se traduce). Recuadros = codones; línea punteada = codón partido entre exones (se completa con el exón vecino).") +
    tema_base(12) + theme(axis.text = element_blank())
  mostrar(p7)
  guardar_figura(p7, "fig7_relacion_adn_arn_proteina", 14, 1.4 + 1.8 * length(filas_u))

  if (!is.null(bloques) && nrow(bloques)) {
    cat("  Codones del exón", ex$k, "y sus aminoácidos:\n  ")
    cat(paste0(gsub("T", "U", bloques$triplete), "=", aa3(bloques$aa), ifelse(bloques$completo, "", "*")), sep = " ")
    cat("\n  (* = codón partido entre exones)\n")
    fase <- (idx[!is.na(idx)][1] - 1) %% 3
    cat("  El exón empieza en la posición", idx[!is.na(idx)][1], "de la CDS -> fase", fase,
        "(0 = el exón empieza en la 1.ª base de un codón; 1 = en la 2.ª; 2 = en la 3.ª).\n\n")
  } else cat("  Este exón no tiene región codificante (es un exón de UTR).\n\n")
  preguntas(
    "¿Cuántos codones completos tiene el exón? ¿Hay codones partidos entre exones?",
    "Si se elimina un nucleótido al inicio del exón, ¿qué pasa con todos los codones siguientes?",
    "Explica con esta figura el concepto de colinealidad gen-proteína (Yanofsky, 1964).",
    "¿Cuántos ATG hay en el 5′UTR antes del codón de inicio? ¿Por qué el ribosoma no empieza en ellos?")
  preguntas_busqueda(
    "Ejecuta uso_codones(datos): ¿qué posición del codón (GC1, GC2 o GC3) cambia más entre especies? Relaciónalo con la degeneración del código.",
    "Abre resultados/*_uso_codones.csv: para la leucina (6 codones), ¿qué codón prefiere cada especie?",
    "Busca la secuencia de Kozak (consenso GCCRCCATGG) con buscar_secuencia(datos, 'RCCATGG'): ¿rodea al ATG de inicio?")
  invisible(datos)
}

# ---------- MÓDULO 5: comparación en un tramo variable -----------------------
GRUPOS_AA <- c(A = "Hidrofóbico", V = "Hidrofóbico", L = "Hidrofóbico", I = "Hidrofóbico",
               M = "Hidrofóbico", F = "Hidrofóbico", W = "Hidrofóbico", P = "Especial (G, P)",
               G = "Especial (G, P)", S = "Polar sin carga", T = "Polar sin carga",
               C = "Polar sin carga", Y = "Polar sin carga", N = "Polar sin carga",
               Q = "Polar sin carga", K = "Básico (+)", R = "Básico (+)", H = "Básico (+)",
               D = "Ácido (-)", E = "Ácido (-)")
COL_GRUPOS <- c("Hidrofóbico" = "#E8A33D", "Polar sin carga" = "#5BA3C9",
                "Básico (+)" = "#2E66B5", "Ácido (-)" = "#C8473A", "Especial (G, P)" = "#9A9A9A")

modulo5_comparacion <- function(datos, ventana = NULL, n_codones = 20) {
  if (is.null(datos$aln)) stop("Primero ejecuta el MÓDULO 1 (alineamiento).")
  caja("MÓDULO 5 · Comparación ADN - ARNm - proteína en un tramo variable",
       "Se toma un tramo de la región codificante donde las especies difieren.",
       "Cambio SINÓNIMO (S): cambia el codón pero NO el aminoácido.",
       "Cambio NO SINÓNIMO (N): cambia el aminoácido.")
  esp <- orden_especies(datos); ref <- datos$ref
  P <- matriz_aln(datos$aln$prot[esp])
  C <- matriz_aln(datos$aln$cds[esp])
  ncol_p <- ncol(P)
  # codón de cada columna proteica
  cod_de <- function(e, j) paste(C[e, (3 * j - 2):(3 * j)], collapse = "")
  difs <- vapply(seq_len(ncol_p), function(j)
    sum(vapply(setdiff(esp, ref), function(e) cod_de(e, j) != cod_de(ref, j), TRUE)), 1)
  pos_ref <- columnas_a_posicion(P[ref, ])
  n_codones <- min(n_codones, ncol_p)
  if (is.null(ventana)) {
    puntaje <- as.numeric(stats::filter(difs, rep(1, n_codones), sides = 1))
    fin <- which.max(puntaje); c1 <- max(1, fin - n_codones + 1)
  } else {
    c1 <- which(pos_ref == ventana[1])[1] %||% 1
    n_codones <- min(25, ventana[2] - ventana[1] + 1)
  }
  c2 <- min(ncol_p, c1 + n_codones - 1)
  cols <- c1:c2
  r1 <- pos_ref[cols][!is.na(pos_ref[cols])]
  cat(sprintf("  Tramo analizado: codones %d a %d de la proteína de %s.\n\n", min(r1), max(r1), nombre_esp(ref)))
  datos$ventana <- c(min(r1), max(r1))

  tiles <- list(); marcas <- list(); conteo <- list()
  y0 <- 0
  for (e in esp) {
    S <- 0; N <- 0
    for (jj in seq_along(cols)) {
      j <- cols[jj]
      cod <- C[e, (3 * j - 2):(3 * j)]; codr <- C[ref, (3 * j - 2):(3 * j)]
      aa <- P[e, j]; aar <- P[ref, j]
      x <- (jj - 1) * 3 + 1:3
      dif_nt <- cod != codr
      tiles[[length(tiles) + 1]] <- data.frame(
        x = x, y = y0, l = cod,
        cat = if (e == ref) cod else ifelse(dif_nt, cod, "igual"), fila = "ADN")
      tiles[[length(tiles) + 1]] <- data.frame(
        x = x, y = y0 - 1, l = ifelse(cod == "T", "U", cod),
        cat = if (e == ref) ifelse(cod == "T", "U", cod) else ifelse(dif_nt, ifelse(cod == "T", "U", cod), "igual"), fila = "ARNm")
      tiles[[length(tiles) + 1]] <- data.frame(
        x = mean(x), y = y0 - 2.05, l = aa3(aa),
        cat = if (e == ref) "aa_ref" else if (aa != aar) "aa_dif" else "aa_igual", fila = "Proteína")
      if (e != ref && any(dif_nt)) {
        tipo <- if (aa == "-" || aar == "-") "Indel" else if (aa == aar) "S" else "N"
        if (tipo == "S") S <- S + 1 else if (tipo == "N") N <- N + 1
        marcas[[length(marcas) + 1]] <- data.frame(x = mean(x), y = y0 + 0.75, tipo = tipo)
      }
    }
    conteo[[e]] <- data.frame(Especie = nombre_esp(e), S = S, N = N)
    y0 <- y0 - 4.2
  }
  dt <- do.call(rbind, tiles)
  dm <- if (length(marcas)) do.call(rbind, marcas) else NULL
  ws <- ifelse(dt$fila == "Proteína", 2.85, 0.92)
  dt$w <- ws
  et <- data.frame(y = unlist(lapply(seq_along(esp), function(i) c(0, -1, -2.05) - (i - 1) * 4.2)),
                   lab = rep(c("ADN", "ARNm", "Proteína"), length(esp)))
  esp_lab <- data.frame(y = 0.75 - (seq_along(esp) - 1) * 4.2, lab = nombre_esp(esp))
  fills <- c(COL_NT, igual = "#EDEDED", aa_ref = "#F6C343", aa_igual = "#F3EBD3", aa_dif = "#E8833A")
  p8 <- ggplot(dt) +
    geom_tile(aes(x, y, fill = cat, width = w), height = 0.78, colour = "white") +
    geom_text(aes(x, y, label = l, colour = cat %in% c("igual", "aa_igual", "aa_ref")),
              size = ifelse(dt$fila == "Proteína", 3, 3.3), fontface = "bold") +
    geom_text(data = et, aes(x = -0.8, y = y, label = lab), hjust = 1, size = 3.2) +
    geom_text(data = esp_lab, aes(x = -0.8, y = y, label = lab), hjust = 1, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = fills, guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "#1B2A33", `FALSE` = "white"), guide = "none") +
    scale_x_continuous(limits = c(-9, length(cols) * 3 + 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 8. Tramo variable de ", datos$gen, ": codones ", min(r1), "-", max(r1), " (numeración de ", nombre_esp(ref), ")"),
         subtitle = paste0("Coloreado = diferente de ", nombre_esp(ref), ". Proteína naranja = aminoácido distinto. ",
                           "Marcas: S = cambio sinónimo, N = no sinónimo, Indel = inserción o deleción.")) +
    tema_base(12) + theme(axis.text = element_blank())
  if (!is.null(dm))
    p8 <- p8 + geom_label(data = dm, aes(x, y, label = tipo, fill = NULL, colour = NULL),
                          size = 2.8, label.padding = unit(0.1, "lines"),
                          colour = ifelse(dm$tipo == "N", "#C04A10", ifelse(dm$tipo == "S", "#1F4E79", "#555555")))
  mostrar(p8)
  guardar_figura(p8, "fig8_tramo_variable", 15, 1.4 + 1.55 * length(esp))

  # Conteo en toda la CDS
  total <- do.call(rbind, lapply(setdiff(esp, ref), function(e) {
    S <- N <- I <- 0
    for (j in seq_len(ncol_p)) {
      if (cod_de(e, j) == cod_de(ref, j)) next
      if (P[e, j] == "-" || P[ref, j] == "-") I <- I + 1 else if (P[e, j] == P[ref, j]) S <- S + 1 else N <- N + 1
    }
    data.frame(Especie = nombre_esp(e), `Codones sinónimos` = S, `Codones no sinónimos` = N,
               `Codones con indel` = I, check.names = FALSE)
  }))
  cat("  Diferencias con", nombre_esp(ref), "en TODA la región codificante:\n")
  print(total, row.names = FALSE)
  datos$tabla_cambios <- total
  write.csv(total, ruta("resultados", paste0(datos$gen, "_cambios_sinonimos.csv")), row.names = FALSE)

  # Figura 8b: composición de aminoácidos (tramo y proteína completa)
  comp <- do.call(rbind, lapply(esp, function(e) {
    tramo <- P[e, cols]; tramo <- tramo[tramo != "-"]
    prot <- letras(datos$modelos[[e]]$proteina)
    rbind(data.frame(Especie = nombre_esp(e), Grupo = GRUPOS_AA[tramo], ambito = "Tramo analizado"),
          data.frame(Especie = nombre_esp(e), Grupo = GRUPOS_AA[prot], ambito = "Proteína completa"))
  }))
  comp <- comp[!is.na(comp$Grupo), ]
  comp$Especie <- factor(comp$Especie, levels = rev(nombre_esp(esp)))
  comp$Grupo <- factor(comp$Grupo, levels = names(COL_GRUPOS))
  comp$ambito <- factor(comp$ambito, levels = c("Tramo analizado", "Proteína completa"))
  p8b <- ggplot(comp, aes(y = Especie, fill = Grupo)) +
    geom_bar(position = "fill", width = 0.7) +
    facet_wrap(~ambito) +
    scale_x_continuous(labels = function(v) paste0(v * 100, "%")) +
    scale_fill_manual(values = COL_GRUPOS, name = "Tipo de aminoácido") +
    labs(x = "Proporción de aminoácidos", y = NULL,
         title = "Figura 8b. Composición de aminoácidos por propiedades químicas",
         subtitle = "¿Los cambios entre especies alteran la composición química o se reemplazan aminoácidos parecidos?") +
    tema_base(12)
  mostrar(p8b)
  guardar_figura(p8b, "fig8b_composicion_aminoacidos", 12, 2 + 0.6 * length(esp))
  tabla_aa <- do.call(rbind, lapply(esp, function(e) {
    t <- table(factor(letras(datos$modelos[[e]]$proteina), levels = names(GRUPOS_AA)))
    data.frame(Especie = nombre_esp(e), t(as.matrix(t)), check.names = FALSE)
  }))
  write.csv(tabla_aa, ruta("resultados", paste0(datos$gen, "_composicion_aminoacidos.csv")), row.names = FALSE)
  preguntas(
    "En el tramo, ¿hay más cambios sinónimos o no sinónimos? ¿Por qué esperarías eso?",
    "¿En qué posición del codón (1.ª, 2.ª o 3.ª) ocurren la mayoría de los cambios sinónimos?",
    "Cuando cambia un aminoácido, ¿se reemplaza por otro del mismo tipo químico (Figura 8b)?",
    "Calcula la razón (no sinónimos / sinónimos) para cada especie. ¿Qué indica sobre la selección en este gen?")
  preguntas_busqueda(
    paste0("En UniProt (", url_uniprot(datos$gen), ") identifica los dominios o sitios activos de la proteína: ",
           "¿los cambios N de la Figura 8c caen fuera de ellos?"),
    paste0("En Ensembl (", url_ensembl(datos$gen), ") abre el 'Gene tree': ¿el árbol del gen coincide con el de las especies?"))
  mapa_cambios(datos)
  invisible(datos)
}

# ---------- MÓDULO 6: efecto de una mutación ---------------------------------
# Tramo más largo de codones CAG/CAA (poliglutamina), p. ej. en HTT
tramo_cag <- function(cds) {
  cod <- substring(cds, seq(1, nchar(cds) - 2, 3), seq(3, nchar(cds), 3))
  r <- rle(cod %in% c("CAG", "CAA")); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  k <- which(r$values)[which.max(r$lengths[r$values])]
  list(ini = ini[k], fin = fin[k], n = r$lengths[k])
}

mutacion_automatica <- function(tipo, cds, n_rep = 25) {
  if (tipo == "auto_cag") {
    t <- tramo_cag(cds)
    if (!length(t$ini)) stop("No se encontró un tramo de repeticiones CAG en la región codificante.")
    return(list(tipo = "insercion", posicion = 3 * t$fin, nuevo = strrep("CAG", n_rep),
                etiqueta = sprintf("El tramo normal tiene %d codones CAG/CAA (codones %d-%d); se añaden %d -> %d repeticiones",
                                   t$n, t$ini, t$fin, n_rep, t$n + n_rep)))
  }
  n <- floor(nchar(cds) / 3)
  if (tipo == "auto_cambio_marco")
    return(list(tipo = "delecion", posicion = 28, longitud = 1,
                etiqueta = "Deleción de 1 nucleótido en el codón 10"))
  for (c in 2:(n - 1)) {
    trip <- substr(cds, 3 * c - 2, 3 * c)
    aa <- traducir(trip)
    for (p in if (tipo == "auto_silenciosa") 3 else 1:3) for (b in c("A", "C", "G", "T")) {
      t2 <- trip; substr(t2, p, p) <- b
      if (t2 == trip) next
      a2 <- traducir(t2)
      ok <- if (tipo == "auto_silenciosa") a2 == aa else a2 == "*"
      if (ok) return(list(tipo = "sustitucion", posicion = 3 * c - 3 + p,
                          ref = substr(trip, p, p), nuevo = b,
                          etiqueta = paste0("Codón ", c, ": ", trip, " -> ", t2)))
    }
  }
  stop("No se encontró una mutación automática de ese tipo.")
}

modulo6_mutacion <- function(datos, ejemplo = 1, propia = NULL) {
  gen <- datos$gen
  menu <- menu_mutaciones(gen)
  if (ejemplo == 0) {
    if (is.null(propia)) stop("Define MUTACION_PROPIA o elige un número de MUTACION_EJEMPLO.")
    if (!is.null(propia$hgvs) && nzchar(propia$hgvs)) {
      mut <- parsear_hgvs(propia$hgvs); mut$etiqueta <- paste("Mutación propia:", propia$hgvs)
    } else mut <- list(tipo = propia$tipo, posicion = propia$posicion_cds, nuevo = toupper(propia$nuevo %||% ""),
                       longitud = propia$longitud %||% 1, etiqueta = "Mutación propia")
    especie <- resolver_especies(propia$especie %||% datos$ref)
  } else {
    if (ejemplo > length(menu)) stop("MUTACION_EJEMPLO debe ser un número entre 0 y ", length(menu))
    mut <- menu[[ejemplo]]
    especie <- mut$especie %||% datos$ref
    if (!is.null(mut$hgvs) && is.null(mut$posicion)) { et <- mut$etiqueta; mut <- parsear_hgvs(mut$hgvs); mut$etiqueta <- et }
  }
  m <- datos$modelos[[especie]]
  if (is.null(m)) stop("La especie elegida para la mutación no está entre las cargadas.")
  cds <- m$cds_sec
  if (startsWith(mut$tipo, "auto_")) {
    et <- mut$etiqueta
    mut <- mutacion_automatica(mut$tipo, cds, mut$n %||% 25)
    mut$etiqueta <- paste0(et, " - ", mut$etiqueta)
  }
  caja(paste0("MÓDULO 6 · Efecto de una mutación en ", gen, " (", nombre_esp(especie), ")"),
       mut$etiqueta,
       "Se aplica el cambio en la región codificante (posición 1 = A del ATG) y se",
       "comparan el ADN, el ARNm y la proteína normales y mutantes.")
  # secuencia desde el ATG hasta después del final del transcrito (para cambios de marco)
  extra <- substr(m$sec, m$tes + 1, min(m$L, m$tes + 300))
  normal <- paste0(substr(m$arnm, sum(m$mapa_arn < m$atg) + 1, nchar(m$arnm)), extra)
  p <- mut$posicion
  if (p < 1 || p > nchar(cds)) stop("La posición debe estar entre 1 y ", nchar(cds), " (longitud de la CDS).")
  aviso <- NULL
  revisar_ref <- function(a, b, esperado) {
    obs <- substr(normal, a, b)
    if (!is.null(esperado) && nzchar(esperado) && esperado != obs)
      aviso <<- sprintf("Atención: en esta secuencia c.%d%s es %s y no %s (puede ser otro transcrito o especie).",
                        a, if (b > a) paste0("_", b) else "", obs, esperado)
    obs
  }
  largo <- mut$longitud %||% 1
  s_ini <- p; s_fin <- p
  if (mut$tipo == "sustitucion") {
    ref_obs <- revisar_ref(p, p, mut$ref)
    mutante <- normal; substr(mutante, p, p) <- mut$nuevo
    cambio_txt <- sprintf("c.%d%s>%s", p, ref_obs, mut$nuevo)
    delta <- 0
  } else if (mut$tipo == "delecion") {
    obs <- revisar_ref(p, p + largo - 1, mut$ref)
    mutante <- paste0(substr(normal, 1, p - 1), substr(normal, p + largo, nchar(normal)))
    cambio_txt <- if (largo == 1) sprintf("c.%ddel%s", p, obs) else sprintf("c.%d_%ddel%s", p, p + largo - 1, if (largo <= 6) obs else "")
    delta <- -largo; s_fin <- p + largo - 1
  } else if (mut$tipo == "duplicacion") {
    dup <- revisar_ref(p, p + largo - 1, mut$ref)
    mutante <- paste0(substr(normal, 1, p + largo - 1), dup, substr(normal, p + largo, nchar(normal)))
    cambio_txt <- if (largo == 1) sprintf("c.%ddup%s", p, dup) else sprintf("c.%d_%ddup%s", p, p + largo - 1, if (largo <= 6) dup else "")
    delta <- largo; s_fin <- p + largo - 1
  } else if (mut$tipo == "insercion") {
    ins <- toupper(mut$nuevo)
    mutante <- paste0(substr(normal, 1, p), ins, substr(normal, p + 1, nchar(normal)))
    cambio_txt <- sprintf("c.%d_%dins%s", p, p + 1, ins)
    delta <- nchar(ins); s_fin <- p + 1
  } else if (mut$tipo == "delins") {
    ins <- toupper(mut$nuevo)
    mutante <- paste0(substr(normal, 1, p - 1), ins, substr(normal, p + largo, nchar(normal)))
    cambio_txt <- sprintf("c.%d%sdelins%s", p, if (largo > 1) paste0("_", p + largo - 1) else "", ins)
    delta <- nchar(ins) - largo; s_fin <- p + largo - 1
  } else stop("Tipo de mutación no reconocido: usa 'sustitucion', 'delecion', 'insercion' o la notación HGVS (hgvs = 'c.68_69del').")
  if (!is.null(aviso)) cat("  ", aviso, "\n\n")

  pn <- traducir_hasta_stop(normal); pm <- traducir_hasta_stop(mutante)
  an <- letras(pn); am <- letras(pm)
  k <- which(an[seq_len(min(length(an), length(am)))] != am[seq_len(min(length(an), length(am)))])[1]
  if (is.na(k) && length(an) != length(am)) k <- min(length(an), length(am)) + 1
  codon_mut <- (p - 1) %/% 3 + 1
  if (is.na(k)) {
    clase <- "Silenciosa (sinónima): el aminoácido no cambia"
    prot_txt <- sprintf("p.%s%d=", aa3(an[codon_mut]), codon_mut)
  } else if (k == 1) {
    clase <- "Pérdida del codón de inicio: no se inicia la traducción normal"
    prot_txt <- "p.Met1?"
  } else if (delta %% 3 != 0) {
    clase <- "Cambio del marco de lectura: todos los codones siguientes cambian"
    prot_txt <- sprintf("p.%s%d%sfsTer%d", aa3(an[k]), k, aa3(am[k]), length(am) - k + 1)
  } else if (am[k] == "*") {
    clase <- "Sin sentido: aparece un codón de parada prematuro (proteína truncada)"
    prot_txt <- sprintf("p.%s%dTer", aa3(an[k]), k)
  } else if (length(am) > length(an) && k == length(an)) {
    clase <- "Pérdida del codón de parada: la proteína se alarga"
    prot_txt <- sprintf("p.Ter%d%sext", k, aa3(am[k]))
  } else if (delta != 0) {
    clase <- sprintf("%s en marco: se %s %d aminoácido(s)", if (delta < 0) "Deleción" else "Inserción",
                     if (delta < 0) "pierde(n)" else "añade(n)", abs(delta) / 3)
    d <- abs(delta) / 3
    prot_txt <- if (delta < 0)
      sprintf("p.%s%d_%s%ddelins%s", aa3(an[k]), k, aa3(an[k + d]), k + d, paste(aa3(am[k]), collapse = ""))
    else sprintf("p.%s%ddelins%s", aa3(an[k]), k, paste(aa3(am[k:(k + d)]), collapse = ""))
  } else {
    clase <- "De sentido erróneo (missense): cambia un aminoácido"
    prot_txt <- sprintf("p.%s%d%s", aa3(an[k]), k, aa3(am[k]))
  }
  cat("  Cambio en el ADN:      ", cambio_txt, "\n")
  cat("  Cambio en la proteína: ", prot_txt, "\n")
  cat("  Clasificación:         ", clase, "\n")
  largo_mut <- if (!is.na(k) && k == 1) 0 else sum(am != "*")
  # Degradación del ARNm mediada por mutaciones sin sentido (NMD), regla de los 50-55 nt
  nmd <- NA
  ptc <- !is.na(k) && k > 1 && am[length(am)] == "*" &&
    ((delta %% 3 != 0 && length(am) < length(an)) || am[k] == "*")
  if (ptc) {
    utr5 <- sum(m$mapa_arn < m$atg)
    fin_ex <- cumsum(m$exones$fin - m$exones$inicio + 1) - utr5
    fin_ex <- fin_ex[-length(fin_ex)]
    if (!length(fin_ex) || max(fin_ex) <= 0) {
      nmd <- "Sin uniones exón-exón después del codón de parada: el ARNm escapa al NMD (posible proteína truncada)."
    } else {
      ultima <- max(fin_ex); if (ultima >= p) ultima <- ultima + delta
      pos_ptc <- 3 * length(am) - 2
      dist <- ultima - pos_ptc
      nmd <- if (dist > 50) sprintf("PROBABLE NMD: el codón de parada prematuro está %d nt antes de la última unión exón-exón (> 50-55 nt): el ARNm se degradaría y casi no habría proteína.", dist)
             else if (dist <= 0) "ESCAPA al NMD: el codón de parada prematuro está en el último exón: podría producirse una proteína truncada."
             else sprintf("ESCAPA al NMD: el codón de parada prematuro está a solo %d nt de la última unión exón-exón (< 50 nt): podría producirse una proteína truncada.", dist)
    }
    cat("  NMD:", nmd, "\n")
  }
  cat("  Longitud de la proteína: normal =", sum(an != "*"), "aa | mutante =",
      if (largo_mut == 0) "0 aa (sin codón de inicio no se produce la proteína normal)" else paste(largo_mut, "aa"), "\n")
  if (gen == "HBB" && !is.na(k))
    cat("  (En la β-globina la numeración clásica no cuenta la Met inicial: codón", k, "=", k - 1, "clásico.)\n")
  cat("\n")
  datos$mutacion <- list(etiqueta = mut$etiqueta, especie = especie, adn = cambio_txt,
                         proteina = prot_txt, clase = clase,
                         largo_normal = sum(an != "*"), largo_mutante = largo_mut,
                         codon = if (is.na(k)) codon_mut else k, nmd = nmd, delta = delta, posicion = p,
                         ptc = if (ptc) length(am) else NA,
                         etiqueta_aa = if (is.na(k)) aa3(an[codon_mut]) else paste0(aa3(an[k]), k, if (k <= length(am)) aa3(am[k]) else ""))

  # Figura 9
  c1 <- max(1, codon_mut - 6); c2 <- c1 + 14
  rango <- (3 * c1 - 2):(3 * c2)
  ln <- letras(substr(normal, min(rango), max(rango)))
  lm <- letras(substr(mutante, min(rango), max(rango)))
  lm <- c(lm, rep(" ", length(ln) - length(lm)))[seq_along(ln)]
  xs <- seq_along(ln)
  cambio_nt <- ln != lm
  filas <- list(
    data.frame(x = xs, y = 6, l = ln, cat = ifelse(xs + min(rango) - 1 >= s_ini & xs + min(rango) - 1 <= s_fin, "sitio", ln)),
    data.frame(x = xs, y = 5, l = lm, cat = ifelse(cambio_nt, "cambio", lm)),
    data.frame(x = xs, y = 3.7, l = gsub("T", "U", ln), cat = gsub("T", "U", ln)),
    data.frame(x = xs, y = 2.7, l = gsub("T", "U", lm), cat = ifelse(cambio_nt, "cambio", gsub("T", "U", lm))))
  dn <- do.call(rbind, filas)
  dn$cat[dn$l == " "] <- "vacio"
  ncod <- c2 - c1 + 1
  aa_n <- an[c1:c2]; aa_m <- am[c1:min(c2, length(am))]
  aa_m <- c(aa_m, rep("", ncod - length(aa_m)))
  dp <- rbind(
    data.frame(x = (seq_len(ncod) - 1) * 3 + 2, y = 1.4, l = aa3(aa_n), cat = ifelse(aa_n == "*", "stop", "aa")),
    data.frame(x = (seq_len(ncod) - 1) * 3 + 2, y = 0.4, l = ifelse(aa_m == "", "", aa3(aa_m)),
               cat = ifelse(aa_m == "*", "stop", ifelse(aa_m == "", "vacio", ifelse(aa_m != aa_n, "aa_cambio", "aa")))))
  et <- data.frame(y = c(6, 5, 3.7, 2.7, 1.4, 0.4),
                   lab = c("ADN normal", "ADN mutante", "ARNm normal", "ARNm mutante", "Proteína normal", "Proteína mutante"))
  num <- data.frame(x = (seq_len(ncod) - 1) * 3 + 2, y = 6.7, lab = c1:c2)
  fills <- c(COL_NT, sitio = "#7B3294", cambio = "#7B3294", aa = "#F3EBD3", aa_cambio = "#E8833A", stop = "#D62728", vacio = "#FFFFFF")
  p9a <- ggplot() +
    geom_tile(data = dn, aes(x, y, fill = cat), width = 0.92, height = 0.7, colour = "white") +
    geom_text(data = dn, aes(x, y, label = l), colour = "white", fontface = "bold", size = 3.3) +
    geom_tile(data = dp, aes(x, y, fill = cat), width = 2.85, height = 0.7, colour = "white") +
    geom_text(data = dp, aes(x, y, label = l, colour = cat == "aa"), fontface = "bold", size = 3.2) +
    geom_text(data = et, aes(x = -0.5, y = y, label = lab), hjust = 1, size = 3.4) +
    geom_text(data = num, aes(x, y, label = lab), size = 2.8, colour = "#777777") +
    geom_text(aes(x = -0.5, y = 6.7, label = "Codón n.º"), hjust = 1, size = 2.8, colour = "#777777") +
    scale_fill_manual(values = fills, guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "#1B2A33", `FALSE` = "white"), guide = "none") +
    scale_x_continuous(limits = c(-10, length(ln) + 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + tema_base(12) + theme(axis.text = element_blank())
  barras <- data.frame(prot = factor(c("Normal", "Mutante"), levels = c("Mutante", "Normal")),
                       aa = c(sum(an != "*"), largo_mut))
  p9b <- ggplot(barras, aes(aa, prot, fill = prot)) + geom_col(width = 0.6) +
    geom_text(aes(label = paste(aa, "aa")), hjust = -0.15, size = 3.6) +
    scale_fill_manual(values = c(Normal = "#1F4E79", Mutante = "#E8833A"), guide = "none") +
    scale_x_continuous(limits = c(0, max(barras$aa) * 1.15)) +
    labs(x = "Longitud de la proteína (aminoácidos)", y = NULL) + tema_base(12)
  p9 <- (p9a / p9b) + plot_layout(heights = c(4, 1)) +
    plot_annotation(title = paste0("Figura 9. Efecto de la mutación ", cambio_txt, " en ", gen, " (", nombre_esp(especie), ")"),
                    subtitle = paste0(prot_txt, " · ", clase, if (!is.na(nmd)) paste0("\n", nmd) else ""),
                    theme = theme(plot.title = element_text(face = "bold", size = 15)))
  mostrar(p9)
  guardar_figura(p9, "fig9_efecto_mutacion", 14, 6.5)
  preguntas(
    "¿Qué tipo de mutación es y por qué? Usa la Figura 9 para justificarlo.",
    "Prueba otra mutación del menú (o una propia). ¿Cuál tiene el efecto más grave en la proteína?",
    "¿Por qué una deleción de 3 nt suele ser menos grave que una de 1 o 2 nt?",
    "Si aparece un codón de parada prematuro, ¿el script predice NMD? Explica por qué la posición del codón importa.",
    "¿El aminoácido afectado es el mismo en las otras especies (Figura 8 o acercamiento(datos, 'c.N'))? ¿Qué sugiere eso?")
  preguntas_busqueda(
    paste0("Busca la variante en ClinVar (", url_clinvar(gen), ") escribiendo su nombre (p. ej. ", cambio_txt,
           "): ¿cuál es su clasificación clínica (patogénica, benigna, incierta)?"),
    paste0("En OMIM (", url_omim(gen), ") lee la descripción de la enfermedad asociada: ¿es dominante o recesiva?"))
  invisible(datos)
}


# -----------------------------------------------------------------------------
# 6. REVISIÓN DE RESPUESTAS
# -----------------------------------------------------------------------------

CLAVES <- list(
  codon_inicio = list(validas = "ATG", nombre = "Codón de inicio",
    pista = "Codifica metionina y marca el inicio del marco de lectura. En el ARNm se escribe AUG."),
  codones_stop = list(validas = c("TAA", "TAG", "TGA"), conjunto = TRUE, nombre = "Codones de parada",
    pista = "Son tres y no codifican ningún aminoácido. En el ARNm: UAA, UAG y UGA."),
  caja_TATA = list(validas = c("TATAAA", "TATAAAA", "TATAAAT", "TATATAA", "TATATA", "TATAWAW", "TATAWA", "TATAWAWR", "TATAAA"),
    nombre = "Caja TATA",
    pista = "Secuencia rica en T y A, unos 25-30 nt antes del TSS. Empieza por TATA... (consenso de 6-7 nt)."),
  caja_CCAAT = list(validas = c("CCAAT", "CCAAAT", "GGCCAATCT", "CCAAT"), nombre = "Caja CCAAT",
    pista = "Su nombre es su secuencia. Se encuentra unos 60-100 nt antes del TSS."),
  sitio_donador = list(validas = c("GT", "GTAAGT", "GTRAGT", "GTAAG", "GTGAGT", "GTAGGT", "GTAAGA", "GTRAG"),
    nombre = "Sitio donador (inicio del intrón)",
    pista = "Regla GT-AG: escribe el dinucleótido con el que empiezan casi todos los intrones."),
  sitio_aceptor = list(validas = c("AG", "CAG", "TAG", "YAG", "NCAG", "YAG"), nombre = "Sitio aceptor (final del intrón)",
    pista = "Regla GT-AG: escribe el dinucleótido con el que terminan casi todos los intrones."),
  senal_poliA = list(validas = c("AATAAA", "ATTAAA"), nombre = "Señal de poliadenilación",
    pista = "Hexámero rico en A, 10-30 nt antes del final del ARNm, donde se añade la cola poli(A).")
)

normalizar <- function(x) {
  x <- toupper(gsub("[^A-Za-z,;]", "", paste(x, collapse = ",")))
  x <- gsub("U", "T", x)
  v <- unlist(strsplit(x, "[,;]"))
  v[nzchar(v)]
}

buscar_patron <- function(s, patron, max_dif = 0) {
  h <- matchPattern(patron, DNAString(s), fixed = FALSE, max.mismatch = max_dif)
  start(h)
}

revisar_respuestas <- function(datos, respuestas) {
  caja("REVISIÓN DE RESPUESTAS",
       "Por cada respuesta correcta se busca esa secuencia en los genes y se grafica.")
  esp <- orden_especies(datos); mods <- datos$modelos[esp]; nombres <- nombre_esp(esp)
  puntos <- 0; total <- 0; graficos <- list(); registro <- list()
  for (campo in names(CLAVES)) {
    cl <- CLAVES[[campo]]
    total <- total + 1
    r <- normalizar(respuestas[[campo]] %||% "")
    if (!length(r)) {
      cat(sprintf("  ○ %-36s sin responder\n", cl$nombre))
      registro[[campo]] <- c(cl$nombre, "", "sin responder")
      next
    }
    if (isTRUE(cl$conjunto)) {
      aciertos <- intersect(r, cl$validas)
      nota <- length(aciertos) / length(cl$validas)
      ok <- nota == 1 && all(r %in% cl$validas)
    } else {
      ok <- r[1] %in% cl$validas; nota <- as.numeric(ok)
    }
    puntos <- puntos + nota
    if (ok) {
      cat(sprintf("  ✔ %-36s %s  - ¡Correcto!\n", cl$nombre, paste(r, collapse = ", ")))
    } else if (nota > 0) {
      cat(sprintf("  ◐ %-36s %s  - parcialmente correcto. Pista: %s\n", cl$nombre, paste(r, collapse = ", "), cl$pista))
    } else {
      cat(sprintf("  ✘ %-36s %s  - Revisa. Pista: %s\n", cl$nombre, paste(r, collapse = ", "), cl$pista))
    }
    registro[[campo]] <- c(cl$nombre, paste(r, collapse = ", "), if (ok) "correcto" else if (nota > 0) "parcial" else "incorrecto")
    if (nota == 0) next
    patron <- if (isTRUE(cl$conjunto)) intersect(r, cl$validas) else r[1]

    # --- Gráfico de búsqueda para cada tipo de respuesta ---
    if (campo == "codon_inicio") {
      d <- do.call(rbind, lapply(seq_along(esp), function(i) {
        m <- mods[[i]]; h <- buscar_patron(m$arnm, patron)
        if (!length(h)) return(NULL)
        utr5 <- sum(m$mapa_arn < m$atg)
        data.frame(especie = nombres[i], x = h, tipo = ifelse(h == utr5 + 1, "ATG de inicio (anotado)", "Otro ATG"))
      }))
      lon <- data.frame(especie = nombres, x = vapply(mods, function(m) nchar(m$arnm), 1))
      g0 <- scale_y_discrete(limits = rev(nombres))
      g <- ggplot() + geom_segment(data = lon, aes(x = 1, xend = x, y = especie, yend = especie), colour = "#BBBBBB", linewidth = 2) +
        geom_point(data = d, aes(x, especie, colour = tipo, size = tipo), shape = 124) +
        scale_colour_manual(values = c("ATG de inicio (anotado)" = "#2E9E44", "Otro ATG" = "#9E9E9E"), name = NULL) +
        scale_size_manual(values = c("ATG de inicio (anotado)" = 9, "Otro ATG" = 5), name = NULL) +
        g0 + labs(x = "Posición en el ARNm maduro (nt desde el TSS)", y = NULL, title = "Codón de inicio: ATG en los ARNm",
             subtitle = "No todos los ATG inician la traducción: el ribosoma usa el primero en buen contexto, tras el 5′UTR.")
    } else if (campo == "codones_stop") {
      d <- data.frame(especie = nombres, stop = vapply(mods, function(m) {
        n <- nchar(m$cds_sec); substr(m$cds_sec, n - 2, n) }, ""))
      d <- merge(expand.grid(especie = nombres, codon = c("TAA", "TAG", "TGA")), d)
      d$usa <- d$codon == d$stop
      d$buscado <- d$codon %in% patron
      g <- ggplot(d, aes(codon, especie)) +
        geom_tile(aes(fill = usa), colour = "white", linewidth = 1) +
        geom_text(aes(label = ifelse(usa, "✔ usado", ifelse(buscado, "", "(no lo escribiste)"))), size = 3.4) +
        scale_fill_manual(values = c(`TRUE` = "#F6C343", `FALSE` = "#F1F1F1"), guide = "none") +
        scale_y_discrete(limits = rev(nombres)) +
        labs(x = NULL, y = NULL, title = "Codón de parada usado por cada especie",
             subtitle = "¿Usan todas las especies el mismo codón de parada?")
    } else if (campo %in% c("caja_TATA", "caja_CCAAT")) {
      zona <- if (campo == "caja_TATA") c(-35, -20) else c(-110, -55)
      d <- do.call(rbind, lapply(seq_along(esp), function(i) {
        m <- mods[[i]]; a <- max(1, m$tss - 200); s <- substr(m$sec, a, m$tss + 20)
        h <- buscar_patron(s, patron)
        h1 <- if (nchar(patron) >= 6) setdiff(buscar_patron(s, patron, 1), h) else integer()
        out <- NULL
        if (length(h)) out <- data.frame(especie = nombres[i], x = h + a - 1 - m$tss, tipo = "idéntica")
        if (length(h1)) out <- rbind(out, data.frame(especie = nombres[i], x = h1 + a - 1 - m$tss, tipo = "variante (1 diferencia)"))
        if (campo == "caja_CCAAT") {
          h2 <- buscar_patron(s, revcomp(patron))
          if (length(h2)) out <- rbind(out, data.frame(especie = nombres[i], x = h2 + a - 1 - m$tss, tipo = "en la hebra opuesta (ATTGG)"))
        }
        out
      }))
      g <- ggplot() +
        annotate("rect", xmin = zona[1], xmax = zona[2], ymin = -Inf, ymax = Inf, fill = "#DCE9F2") +
        geom_blank(data = data.frame(especie = nombres), aes(y = especie)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        annotate("text", x = 2, y = Inf, label = "TSS", vjust = 1.5, hjust = 0, size = 3.2) +
        scale_y_discrete(limits = rev(nombres)) +
        coord_cartesian(xlim = c(-200, 20)) +
        labs(x = "Posición relativa al TSS (nt)", y = NULL,
             title = paste0(cl$nombre, ": búsqueda de '", patron, "' en el promotor"),
             subtitle = "Franja azul = zona esperada. Las cajas reales varían: también se marcan variantes con 1 diferencia.")
      if (!is.null(d)) g <- g + geom_point(data = d, aes(x, especie, shape = tipo), size = 4, colour = "#C8473A") +
        scale_shape_manual(values = c("idéntica" = 16, "variante (1 diferencia)" = 1, "en la hebra opuesta (ATTGG)" = 17), name = NULL)
      else g <- g + annotate("text", x = -90, y = length(nombres) / 2 + 0.5, label = "No se encontró en ninguna especie", colour = "#777777")
    } else if (campo %in% c("sitio_donador", "sitio_aceptor")) {
      d <- do.call(rbind, lapply(seq_along(esp), function(i) {
        m <- mods[[i]]
        if (!nrow(m$intrones)) return(NULL)
        do.call(rbind, lapply(seq_len(nrow(m$intrones)), function(k) {
          a <- m$intrones$inicio[k]; b <- m$intrones$fin[k]; n <- nchar(patron)
          trozo <- if (campo == "sitio_donador") substr(m$sec, a, a + n - 1) else substr(m$sec, b - n + 1, b)
          ok <- length(buscar_patron(trozo, patron)) > 0
          vista <- if (campo == "sitio_donador") paste0(substr(m$sec, a - 3, a - 1), "|", tolower(substr(m$sec, a, a + 5)))
                   else paste0(tolower(substr(m$sec, b - 5, b)), "|", substr(m$sec, b + 1, b + 3))
          data.frame(especie = nombres[i], intron = paste("Intrón", k), ok = ok, vista = vista)
        }))
      }))
      if (is.null(d)) { cat("    (El gen no tiene intrones: no hay sitios de splicing que buscar.)\n"); next }
      g <- ggplot(d, aes(intron, especie)) + geom_tile(aes(fill = ok), colour = "white", linewidth = 1) +
        geom_text(aes(label = vista), family = "mono", size = 3.6) +
        scale_fill_manual(values = c(`TRUE` = "#BFE3C0", `FALSE` = "#F4C7C3"), labels = c(`TRUE` = "coincide", `FALSE` = "no coincide"), name = NULL) +
        scale_y_discrete(limits = rev(nombres)) +
        labs(x = NULL, y = NULL, title = paste0(cl$nombre, ": '", patron, "' en cada límite"),
             subtitle = "MAYÚSCULAS = exón, minúsculas = intrón, | = límite de splicing.")
    } else if (campo == "senal_poliA") {
      d <- do.call(rbind, lapply(seq_along(esp), function(i) {
        m <- mods[[i]]; a <- max(1, m$tes - 60); s <- substr(m$sec, a, min(m$L, m$tes + 10))
        h <- buscar_patron(s, patron)
        if (!length(h)) return(NULL)
        data.frame(especie = nombres[i], x = h + a - 1 - m$tes)
      }))
      g <- ggplot() + annotate("rect", xmin = -35, xmax = -10, ymin = -Inf, ymax = Inf, fill = "#DCE9F2") +
        geom_blank(data = data.frame(especie = nombres), aes(y = especie)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        annotate("text", x = 1, y = Inf, label = "fin del ARNm", vjust = 1.5, hjust = 0, size = 3.2) +
        scale_y_discrete(limits = rev(nombres)) + coord_cartesian(xlim = c(-60, 10)) +
        labs(x = "Posición relativa al final del transcrito (nt)", y = NULL,
             title = paste0("Señal de poliadenilación: '", patron, "'"),
             subtitle = "Franja azul = zona esperada (10-35 nt antes del sitio de corte).")
      if (!is.null(d)) g <- g + geom_point(data = d, aes(x, especie), size = 4, shape = 18, colour = "#7B3294")
    }
    graficos[[campo]] <- g + tema_base(11)
  }
  cat(sprintf("\n  PUNTAJE: %.1f de %d\n\n", puntos, total))
  datos$respuestas <- list(puntos = puntos, total = total, registro = registro)
  if (length(graficos)) {
    p10 <- wrap_plots(graficos, ncol = 2) +
      plot_annotation(title = "Figura 10. Búsqueda de tus respuestas en las secuencias",
                      theme = theme(plot.title = element_text(face = "bold", size = 15)))
    mostrar(p10)
    guardar_figura(p10, "fig10_busqueda_respuestas", 15, 2 + 3.2 * ceiling(length(graficos) / 2))
  }
  invisible(datos)
}


# -----------------------------------------------------------------------------
# 7. INFORME FINAL
# -----------------------------------------------------------------------------

generar_informe <- function(datos, nombre = "Estudiante") {
  archivo <- ruta("resultados", paste0("informe_", datos$gen, "_", gsub("[^A-Za-z0-9]+", "_", iconv(nombre, to = "ASCII//TRANSLIT")), ".txt"))
  con <- file(archivo, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  w <- function(...) writeLines(paste0(...), con)
  w("INFORME - Actividad: Del gen a la proteína en varias especies")
  w("Estudiante: ", nombre, "    Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  w("Gen: ", datos$gen, " - ", CATALOGO[[datos$gen]]$nombre)
  w(strrep("=", 70))
  w("\n1. SECUENCIAS USADAS (NCBI)")
  capture <- capture.output(print(datos$tabla, row.names = FALSE)); w(paste(capture, collapse = "\n"))
  if (!is.null(datos$resumen_identidad)) {
    w("\n2. IDENTIDAD CON ", toupper(nombre_esp(datos$ref)), " (%)")
    w(paste(capture.output(print(datos$resumen_identidad, row.names = FALSE)), collapse = "\n"))
  }
  if (!is.null(datos$tabla_estructura)) {
    w("\n3. ESTRUCTURA DEL GEN")
    w(paste(capture.output(print(datos$tabla_estructura, row.names = FALSE)), collapse = "\n"))
  }
  if (!is.null(datos$tabla_cambios)) {
    w("\n4. CAMBIOS SINÓNIMOS Y NO SINÓNIMOS (toda la CDS)")
    w("Tramo analizado en la Figura 8: codones ", datos$ventana[1], "-", datos$ventana[2])
    w(paste(capture.output(print(datos$tabla_cambios, row.names = FALSE)), collapse = "\n"))
  }
  if (!is.null(datos$mutacion)) {
    mu <- datos$mutacion
    w("\n5. MUTACIÓN ANALIZADA")
    w(mu$etiqueta, " (", nombre_esp(mu$especie), ")")
    w("ADN: ", mu$adn, "   Proteína: ", mu$proteina)
    w("Clasificación: ", mu$clase)
    w("Longitud: normal ", mu$largo_normal, " aa / mutante ", mu$largo_mutante, " aa")
  }
  if (!is.null(datos$respuestas)) {
    w("\n6. RESPUESTAS")
    for (r in datos$respuestas$registro) w(sprintf("  %-36s %-22s %s", r[1], r[2], r[3]))
    w(sprintf("  PUNTAJE: %.1f de %d", datos$respuestas$puntos, datos$respuestas$total))
  }
  if (!is.null(datos$quiz)) {
    w("\n6b. CUESTIONARIO INTERACTIVO")
    for (r in datos$quiz$registro) w(sprintf("  %-60s %-12s %s", r[1], r[2], r[3]))
    w(sprintf("  PUNTAJE: %d de %d", datos$quiz$puntos, datos$quiz$total))
  }
  if (!is.null(datos$variante)) {
    v <- datos$variante
    w("\n6c. VARIANTE EN CONTEXTO")
    w("Exón: ", v$exon, " · Codón afectado: ", v$codon, " · Conservación del aminoácido: ", v$conservacion)
  }
  w("\n7. FIGURAS GENERADAS (carpeta figuras/)")
  for (f in .estado$figuras) w("  ", f)
  caja("INFORME GENERADO", paste0("resultados/", basename(archivo)),
       "Entrega este archivo junto con las figuras y tus respuestas a las preguntas.")
  invisible(archivo)
}


# -----------------------------------------------------------------------------
# 8. AMPLIACIONES: búsquedas, mapa de cambios (Fig. 8c), exones, HGVS y 3D
# -----------------------------------------------------------------------------

# ---------- Enlaces para las preguntas con búsqueda --------------------------
url_ncbi_gene <- function(gen) paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gen, "%5Bsym%5D+AND+human%5Borgn%5D")
url_clinvar   <- function(gen) paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=", gen, "%5Bgene%5D")
url_uniprot   <- function(gen) paste0("https://www.uniprot.org/uniprotkb?query=gene:", gen, "+AND+organism_id:9606")
url_omim      <- function(gen) paste0("https://www.omim.org/search?search=", gen)
url_alphafold <- function(gen) paste0("https://alphafold.ebi.ac.uk/search/text/", gen)
url_ensembl   <- function(gen) paste0("https://www.ensembl.org/Multi/Search/Results?q=", gen)
url_genoma    <- "https://www.ncbi.nlm.nih.gov/datasets/genome/"

preguntas_busqueda <- function(...) {
  cat("  \U0001F50E Preguntas con búsqueda (respalda tu observación con datos):\n")
  q <- c(...)
  cat(paste0("    B", seq_along(q), ". ", q), sep = "\n")
  cat("\n")
}

# ---------- Búsqueda de una secuencia en todos los genes ----------------------
# patron admite códigos IUPAC (R = A/G, Y = C/T, W = A/T, N = cualquiera...).
buscar_secuencia <- function(datos, patron, max_dif = 0, ambas_hebras = FALSE, promotor = 200) {
  patron <- toupper(gsub("U", "T", gsub("\\s", "", patron)))
  esp <- orden_especies(datos); nombres <- nombre_esp(esp)
  hits <- do.call(rbind, lapply(seq_along(esp), function(i) {
    m <- datos$modelos[[esp[i]]]
    a <- max(1, m$tss - promotor); s <- substr(m$sec, a, m$tes)
    pats <- if (ambas_hebras) unique(c(patron, revcomp(patron))) else patron
    do.call(rbind, lapply(pats, function(pt) {
      h <- start(matchPattern(pt, DNAString(s), fixed = FALSE, max.mismatch = max_dif))
      if (!length(h)) return(NULL)
      g <- h + a - 1
      reg <- m$region[g]
      num <- vapply(g, function(x) {
        if (m$region[x] == "intrón") paste("intrón", which(m$intrones$inicio <= x & m$intrones$fin >= x))
        else if (m$region[x] %in% c("CDS", "5′UTR", "3′UTR")) paste("exón", which(m$exones$inicio <= x & m$exones$fin >= x))
        else "" }, "")
      data.frame(especie = nombres[i], posicion_TSS = g - m$tss, region = reg, donde = num,
                 hebra = if (pt == patron) "+" else "-", secuencia = substring(m$sec, g, g + nchar(pt) - 1))
    }))
  }))
  caja(paste0("Búsqueda de '", patron, "'", if (max_dif) paste0(" (hasta ", max_dif, " diferencia/s)") else ""),
       "Se busca en el promotor (200 nt) y en todo el gen (exones e intrones) de cada especie.")
  if (is.null(hits) || !nrow(hits)) {
    cat("  No se encontró en ninguna especie.\n\n"); return(invisible(NULL))
  }
  tab <- as.data.frame.matrix(table(factor(hits$especie, levels = nombres), hits$region))
  print(tab); cat("\n")
  hits$especie <- factor(hits$especie, levels = rev(nombres))
  lon <- data.frame(especie = factor(nombres, levels = rev(nombres)),
                    x0 = -promotor, x1 = vapply(datos$modelos[esp], function(m) m$tes - m$tss, 1))
  exs <- do.call(rbind, lapply(seq_along(esp), function(i) {
    m <- datos$modelos[[esp[i]]]
    data.frame(especie = factor(nombres[i], levels = rev(nombres)), xmin = m$exones$inicio - m$tss, xmax = m$exones$fin - m$tss)
  }))
  reg_col <- c(COL_REG, "intergénico" = "#999999")
  p <- ggplot() +
    geom_segment(data = lon, aes(x = x0, xend = x1, y = especie, yend = especie), colour = "#BBBBBB", linewidth = 1) +
    geom_rect(data = exs, aes(xmin = xmin, xmax = xmax, ymin = as.numeric(especie) - 0.18, ymax = as.numeric(especie) + 0.18),
              fill = "#8DB8D0") +
    geom_point(data = hits, aes(posicion_TSS, especie, colour = region), size = 3.2, shape = 124, stroke = 2) +
    scale_colour_manual(values = reg_col, name = "Región del hallazgo") +
    labs(x = "Posición relativa al TSS (nt)", y = NULL,
         title = paste0("Búsqueda de '", patron, "' en ", datos$gen),
         subtitle = "Cajas azules = exones. Cada marca es una coincidencia.") + tema_base(12)
  mostrar(p)
  guardar_figura(p, paste0("busqueda_", patron), 13, 2 + 0.6 * length(esp))
  write.csv(hits, ruta("resultados", paste0(datos$gen, "_busqueda_", patron, ".csv")), row.names = FALSE)
  invisible(hits)
}

# ---------- Ver un tramo de la secuencia --------------------------------------
ver_region <- function(datos, desde, hasta, especie = NULL) {
  especie <- especie %||% datos$ref
  m <- datos$modelos[[especie]]
  a <- max(1, m$tss + desde); b <- min(m$L, m$tss + hasta)
  s <- m$nt[a:b]; r <- m$region[a:b]
  s[r %in% c("intrón", "promotor", "intergénico")] <- tolower(s[r %in% c("intrón", "promotor", "intergénico")])
  cat(sprintf("\n  %s %s, posiciones %+d a %+d respecto al TSS (MAYÚSCULAS = exón):\n", datos$gen, nombre_esp(especie), desde, hasta))
  trozos <- split(s, ceiling(seq_along(s) / 60))
  for (k in seq_along(trozos)) cat(sprintf("  %+6d  %s\n", desde + (k - 1) * 60, paste(trozos[[k]], collapse = "")))
  cat("\n")
  invisible(paste(s, collapse = ""))
}

# ---------- Uso de codones y contenido GC -------------------------------------
uso_codones <- function(datos) {
  esp <- orden_especies(datos)
  tab <- do.call(rbind, lapply(esp, function(e) {
    s <- letras(datos$modelos[[e]]$cds_sec)
    n <- floor(length(s) / 3) * 3; s <- s[seq_len(n)]
    pos <- rep(1:3, length.out = n)
    gc <- function(v) round(100 * mean(v %in% c("G", "C")), 1)
    data.frame(Especie = nombre_esp(e), `GC total` = gc(s), GC1 = gc(s[pos == 1]),
               GC2 = gc(s[pos == 2]), GC3 = gc(s[pos == 3]), check.names = FALSE)
  }))
  caja("USO DE CODONES Y CONTENIDO GC DE LA REGIÓN CODIFICANTE",
       "GC1, GC2, GC3 = % de G+C en la 1.ª, 2.ª y 3.ª posición de los codones.",
       "La 3.ª posición suele ser sinónima: refleja sesgos mutacionales y de uso de codones.")
  print(tab, row.names = FALSE); cat("\n")
  largo <- reshape(tab, direction = "long", varying = c("GC total", "GC1", "GC2", "GC3"),
                   v.names = "valor", timevar = "pos", times = c("GC total", "GC1", "GC2", "GC3"), idvar = "Especie")
  largo$Especie <- factor(largo$Especie, levels = nombre_esp(esp))
  p <- ggplot(largo, aes(pos, valor, fill = Especie)) + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "#777777") +
    scale_fill_manual(values = c("#1F4E79", "#7B3294", "#E66101", "#008837", "#2C7BB6", "#A6611A", "#5E3C99")[seq_along(esp)]) +
    labs(x = NULL, y = "% G + C", title = paste0("Contenido GC por posición del codón en ", datos$gen),
         subtitle = "¿Qué posición del codón varía más entre especies? ¿Por qué?") + tema_base(12)
  mostrar(p)
  guardar_figura(p, "uso_codones_GC", 11, 5)
  cods <- do.call(rbind, lapply(esp, function(e) {
    s <- datos$modelos[[e]]$cds_sec; n <- floor(nchar(s) / 3)
    c3 <- substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
    t <- table(c3); data.frame(Especie = nombre_esp(e), Codon = names(t), Aminoacido = aa3(vapply(names(t), traducir, "")), Veces = as.integer(t))
  }))
  write.csv(cods, ruta("resultados", paste0(datos$gen, "_uso_codones.csv")), row.names = FALSE)
  cat("  Tabla completa de codones: resultados/", datos$gen, "_uso_codones.csv\n\n", sep = "")
  invisible(tab)
}

perfil_gc <- function(datos, especie = NULL, ventana = 100, promotor = 200) {
  especie <- especie %||% datos$ref
  m <- datos$modelos[[especie]]
  a <- max(1, m$tss - promotor); pos <- a:m$tes
  s <- m$nt[pos]
  gcv <- as.numeric(s %in% c("G", "C"))
  cg <- c(as.numeric(s[-length(s)] == "C" & s[-1] == "G"), 0)
  cc <- as.numeric(s == "C"); gg <- as.numeric(s == "G")
  f <- function(v) as.numeric(stats::filter(v, rep(1, ventana), sides = 2))
  gc <- 100 * f(gcv) / ventana
  oe <- f(cg) * ventana / pmax(f(cc) * f(gg), 1)
  d <- data.frame(x = pos - m$tss, GC = gc, oe = oe)
  isla <- !is.na(d$GC) & d$GC > 50 & d$oe > 0.6
  reg <- data.frame(x = pos - m$tss, reg = factor(m$region[pos], levels = names(COL_REG)))
  p1 <- ggplot(d, aes(x, GC)) + geom_line(colour = "#1F4E79", na.rm = TRUE) +
    geom_hline(yintercept = 50, linetype = "dashed") + labs(x = NULL, y = "% GC") + tema_base(11)
  p2 <- ggplot(d, aes(x, oe)) + geom_line(colour = "#C8473A", na.rm = TRUE) +
    geom_hline(yintercept = 0.6, linetype = "dashed") + labs(x = NULL, y = "CpG obs/esp") + tema_base(11)
  p3 <- ggplot(reg, aes(x, 1, fill = reg)) + geom_raster() +
    geom_point(data = d[isla, ], aes(x, 1.5), inherit.aes = FALSE, shape = 15, colour = "#7B3294", size = 1.2) +
    scale_fill_manual(values = COL_REG, drop = TRUE, name = "Región") +
    labs(x = "Posición relativa al TSS (nt)", y = NULL) + tema_base(11) + theme(axis.text.y = element_blank())
  p <- (p1 / p2 / p3) + plot_layout(heights = c(2, 2, 1)) +
    plot_annotation(title = paste0("Contenido GC e islas CpG en ", datos$gen, " (", nombre_esp(especie), ")"),
                    subtitle = paste0("Ventana de ", ventana, " nt. Morado = criterio de isla CpG (GC > 50 % y CpG obs/esp > 0,6)."),
                    theme = theme(plot.title = element_text(face = "bold")))
  mostrar(p)
  guardar_figura(p, paste0("perfil_GC_", especie), 13, 6)
  cat("\n  Posiciones con criterio de isla CpG:", sum(isla), "de", sum(!is.na(d$GC)), "\n\n")
  invisible(d)
}

# ---------- Alineamiento de un exón con sus bordes intrónicos -----------------
alinear_exon <- function(datos, exon = 1, flanco = 40, marcar_cds = NULL) {
  esp <- orden_especies(datos); ref <- datos$ref; mr <- datos$modelos[[ref]]
  exon <- max(1, min(exon, nrow(mr$exones)))
  P <- if (!is.null(datos$aln)) matriz_aln(datos$aln$prot[esp]) else NULL
  pos_e <- mr$exones$inicio[exon]:mr$exones$fin[exon]
  idx <- match(pos_e, mr$mapa_cds); idx <- idx[!is.na(idx)]
  piezas <- list(); info <- character()
  for (e in esp) {
    m <- datos$modelos[[e]]
    k <- exon
    if (e != ref && length(idx) && !is.null(P)) {
      cod <- (idx[ceiling(length(idx) / 2)] - 1) %/% 3 + 1
      col <- which(columnas_a_posicion(P[ref, ]) == cod)[1]
      ce <- columnas_a_posicion(P[e, ])[col]
      if (!is.na(ce)) {
        g <- m$mapa_cds[3 * ce - 2]
        k <- which(m$exones$inicio <= g & m$exones$fin >= g)
      }
    }
    if (!length(k) || k > nrow(m$exones)) next
    a <- max(1, m$exones$inicio[k] - flanco); b <- min(m$L, m$exones$fin[k] + flanco)
    s <- m$nt[a:b]
    es_ex <- (a:b) >= m$exones$inicio[k] & (a:b) <= m$exones$fin[k]
    s[!es_ex] <- tolower(s[!es_ex])
    piezas[[e]] <- paste(s, collapse = "")
    info <- c(info, sprintf("%s: exón %d (%d nt)", nombre_esp(e), k, m$exones$fin[k] - m$exones$inicio[k] + 1))
  }
  caja(paste0("Exón ", exon, " de ", datos$gen, " (", nombre_esp(ref), ") y su exón ortólogo en cada especie"), info)
  aln <- alinear(DNAStringSet(toupper(unlist(piezas))))
  # recuperar mayúsculas/minúsculas
  aln <- vapply(names(aln), function(e) {
    orig <- letras(piezas[[e]]); a <- letras(aln[[e]]); k <- 0
    for (j in seq_along(a)) if (a[j] != "-") { k <- k + 1; a[j] <- orig[k] }
    paste(a, collapse = "") }, "")
  M <- do.call(rbind, lapply(aln, letras))
  nombres <- nombre_esp(names(aln))
  ncolm <- ncol(M); ancho <- 100
  d <- data.frame(col = rep(seq_len(ncolm), times = nrow(M)),
                  especie = factor(rep(nombres, each = ncolm), levels = rev(nombres)),
                  l = as.vector(t(M)))
  d$exon <- d$l %in% c("A", "C", "G", "T")
  d$base <- toupper(d$l)
  d$bloque <- paste0("Columnas ", ((d$col - 1) %/% ancho) * ancho + 1, "-", pmin(ncolm, ((d$col - 1) %/% ancho + 1) * ancho))
  d$bloque <- factor(d$bloque, levels = unique(d$bloque))
  d$x <- (d$col - 1) %% ancho + 1
  refc <- M[1, ]
  d$fill <- ifelse(d$base == "-", "-", ifelse(!d$exon, "intrón", d$base))
  p <- ggplot(d, aes(x, especie)) +
    geom_tile(aes(fill = fill, alpha = exon), colour = "white", linewidth = 0.2) +
    geom_text(aes(label = l), size = 2.4, colour = ifelse(d$exon, "white", "#555555")) +
    scale_fill_manual(values = c(COL_NT, "intrón" = "#E9E4D8"), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 1), guide = "none") +
    facet_wrap(~bloque, ncol = 1) +
    labs(x = NULL, y = NULL, title = paste0("Exón ", exon, " de ", datos$gen, " alineado con ", flanco, " nt de intrón a cada lado"),
         subtitle = "MAYÚSCULAS coloreadas = exón; minúsculas sobre beige = intrón; '-' = hueco.") +
    tema_base(10) + theme(axis.text.x = element_blank(), strip.text = element_text(hjust = 0))
  if (!is.null(marcar_cds)) {
    g <- mr$mapa_cds[marcar_cds]
    pr <- g - (max(1, mr$exones$inicio[exon] - flanco)) + 1
    colm <- which(columnas_a_posicion(M[1, ]) == pr)[1]
    if (!is.na(colm)) {
      mk <- data.frame(bloque = levels(d$bloque)[(colm - 1) %/% ancho + 1], x = (colm - 1) %% ancho + 1)
      p <- p + geom_tile(data = mk, aes(x, y = length(nombres) / 2 + 0.5), inherit.aes = FALSE,
                         width = 1, height = length(nombres), fill = NA, colour = "#1B2A33", linewidth = 1)
    }
  }
  mostrar(p)
  guardar_figura(p, paste0("exon", exon, "_alineado"), 15, 1.5 + (0.32 * length(nombres) + 0.5) * nlevels(d$bloque))
  invisible(aln)
}

# ---------- FIGURA 8c: mapa compacto de todos los cambios ----------------------
mapa_cambios <- function(datos, max_codones = 60, por_fila = 22) {
  if (is.null(datos$aln)) stop("Primero ejecuta el MÓDULO 1 (alineamiento).")
  esp <- orden_especies(datos); ref <- datos$ref; mr <- datos$modelos[[ref]]
  P <- matriz_aln(datos$aln$prot[esp]); C <- matriz_aln(datos$aln$cds[esp])
  ncp <- ncol(P); pos_ref <- columnas_a_posicion(P[ref, ])
  cod <- function(e, j) paste(C[e, (3 * j - 2):(3 * j)], collapse = "")
  tipo <- matrix("", length(esp), ncp, dimnames = list(esp, NULL))
  for (e in setdiff(esp, ref)) for (j in seq_len(ncp)) {
    if (cod(e, j) == cod(ref, j)) next
    tipo[e, j] <- if (P[e, j] == "-" || P[ref, j] == "-") "Indel" else if (P[e, j] == "*") "Stop" else if (P[e, j] == P[ref, j]) "S" else "N"
  }
  dif <- which(colSums(tipo != "") > 0)
  omitidos <- 0
  if (length(dif) > max_codones) {
    prio <- vapply(dif, function(j) 10 * sum(tipo[, j] %in% c("N", "Indel", "Stop")) + sum(tipo[, j] == "S"), 1)
    keep <- sort(dif[order(-prio)][seq_len(max_codones)])
    omitidos <- length(dif) - max_codones
    dif <- keep
  }
  stops <- vapply(esp, function(e) { s <- datos$modelos[[e]]$cds_sec; substr(s, nchar(s) - 2, nchar(s)) }, "")
  # exón de cada columna (según la referencia)
  ex_col <- integer(ncp); ult <- 1
  for (j in seq_len(ncp)) {
    if (!is.na(pos_ref[j])) ult <- findInterval(mr$mapa_cds[3 * pos_ref[j] - 2], mr$exones$inicio)
    ex_col[j] <- ult
  }
  utr5 <- sum(mr$mapa_arn < mr$atg); utr3 <- sum(mr$mapa_arn > mr$fin_cds)
  # Construcción de la secuencia de unidades
  U <- list(); x <- 0
  bloque <- function(reg, lab) { U[[length(U) + 1]] <<- list(t = "b", reg = reg, lab = lab, x = x); x <<- x + 2.6 }
  codon <- function(j) { U[[length(U) + 1]] <<- list(t = "c", j = j, x = x); x <<- x + 3.5 }
  bloque("promotor", "Promotor")
  if (utr5 > 0) bloque("5′UTR", paste0("5′UTR ", utr5, " nt"))
  nex <- nrow(mr$exones)
  ex_cds <- unique(ex_col)
  for (k in seq_len(nex)) {
    cols <- which(ex_col == k)
    if (!length(cols)) {
      if (k < min(ex_cds)) { } else if (k > max(ex_cds)) { } else bloque("CDS", paste0("Exón ", k))
    } else {
      corrida <- 0
      for (j in cols) {
        if (j %in% dif) {
          if (corrida > 0) bloque("CDS", paste0(corrida, " cod. iguales"))
          corrida <- 0; codon(j)
        } else corrida <- corrida + 1
      }
      if (corrida > 0) bloque("CDS", paste0(corrida, " cod. iguales"))
    }
    if (k == max(ex_cds)) {
      if (length(unique(stops)) > 1) { U[[length(U) + 1]] <- list(t = "s", x = x); x <- x + 3.5 }
      else bloque("CDS", paste0("Stop ", stops[1]))
    }
    if (k < nex) {
      ln <- mr$intrones$fin[k] - mr$intrones$inicio[k] + 1
      bloque("intrón", paste0("Intrón ", k, " (", ln, " nt)"))
    }
  }
  if (utr3 > 0) bloque("3′UTR", paste0("3′UTR ", utr3, " nt"))

  # Reparte las unidades en filas (tramos) de 'por_fila' codones
  tramo <- 1; ncod <- 0; x <- 0
  for (i in seq_along(U)) {
    if (U[[i]]$t != "b") {
      if (ncod > 0 && ncod %% por_fila == 0) { tramo <- tramo + 1; x <- 0 }
      ncod <- ncod + 1
    }
    U[[i]]$x <- x; U[[i]]$tramo <- tramo
    x <- x + if (U[[i]]$t == "b") 2.6 else 3.5
  }
  xmax <- max(vapply(U, function(u) u$x + if (u$t == "b") 2.6 else 3.5, 1))
  ntr <- tramo
  ne <- length(esp); H <- 2.7
  yc <- function(i) -(i - 1) * H; ya <- function(i) yc(i) - 1
  tiles <- list(); aa_t <- list(); blq <- list(); cab <- list(); exlab <- list()
  for (u in U) {
    if (u$t == "b") {
      blq[[length(blq) + 1]] <- data.frame(x0 = u$x, x1 = u$x + 2.2, reg = u$reg, lab = u$lab, tramo = u$tramo)
    } else {
      for (i in seq_len(ne)) {
        e <- esp[i]
        trip <- if (u$t == "c") C[e, (3 * u$j - 2):(3 * u$j)] else letras(stops[e])
        tref <- if (u$t == "c") C[ref, (3 * u$j - 2):(3 * u$j)] else letras(stops[ref])
        tiles[[length(tiles) + 1]] <- data.frame(x = u$x + 0:2 + 0.5, y = yc(i), l = trip,
          f = if (e == ref) trip else ifelse(trip != tref, trip, "igual"), tramo = u$tramo)
        if (u$t == "c") {
          a <- P[e, u$j]; ty <- tipo[e, u$j]
          f <- if (e == ref) "ref" else if (ty == "") "igual_aa" else paste0("t", ty)
          aa_t[[length(aa_t) + 1]] <- data.frame(x = u$x + 1.5, y = ya(i), l = aa3(a), f = f, tramo = u$tramo)
        } else {
          aa_t[[length(aa_t) + 1]] <- data.frame(x = u$x + 1.5, y = ya(i), l = "Stop",
                                                 f = if (e == ref) "ref" else if (stops[e] == stops[ref]) "igual_aa" else "tStop", tramo = u$tramo)
        }
      }
      cab[[length(cab) + 1]] <- data.frame(x = u$x + 1.5, l = if (u$t == "c") as.character(pos_ref[u$j] %||% "ins") else "fin", tramo = u$tramo)
      if (u$t == "c" && is.na(pos_ref[u$j])) cab[[length(cab)]]$l <- "ins"
    }
  }
  dt <- do.call(rbind, tiles); da <- do.call(rbind, aa_t); db <- do.call(rbind, blq); dc <- do.call(rbind, cab)
  ytop <- yc(1) + 0.6; ybot <- ya(ne) - 0.55
  fills <- c(COL_NT, igual = "#EDEDED", ref = "#F6C343", igual_aa = "#F3EBD3",
             tS = "#9CC3DD", tN = "#E8833A", tIndel = "#9E9E9E", tStop = "#D62728")
  et <- data.frame(y = rep(c(sapply(seq_len(ne), yc), sapply(seq_len(ne), ya)), ntr),
                   lab = rep(c(nombre_esp(esp), rep("aa", ne)), ntr),
                   tramo = rep(seq_len(ntr), each = 2 * ne))
  etiq_tr <- function(v) paste0("Parte ", v, " de ", ntr)
  for (nm in c("db", "dt", "da", "dc", "et")) { z <- get(nm); z$tramo <- factor(etiq_tr(z$tramo), levels = etiq_tr(seq_len(ntr))); assign(nm, z) }
  p <- ggplot() +
    geom_rect(data = db, aes(xmin = x0, xmax = x1, ymin = ybot, ymax = ytop, fill = reg)) +
    geom_text(data = db, aes(x = (x0 + x1) / 2, y = (ytop + ybot) / 2, label = lab,
                             colour = reg %in% c("CDS")), angle = 90, size = 2.9) +
    geom_tile(data = dt, aes(x, y, fill = f), width = 0.95, height = 0.8) +
    geom_text(data = dt, aes(x, y, label = l, colour = f %in% c("igual")), size = 2.7, fontface = "bold") +
    geom_tile(data = da, aes(x, y, fill = f), width = 2.9, height = 0.8) +
    geom_text(data = da, aes(x, y, label = l), size = 2.5) +
    geom_text(data = dc, aes(x, ytop + 0.45, label = l), size = 2.4, colour = "#555555") +
    geom_text(data = et, aes(x = -0.6, y = y, label = lab), hjust = 1,
              size = ifelse(et$lab == "aa", 2.6, 3.4), fontface = ifelse(et$lab == "aa", "plain", "bold")) +
    scale_fill_manual(values = c(fills, COL_REG), breaks = c("tS", "tN", "tIndel", "tStop", "promotor", "5′UTR", "CDS", "intrón"),
                      labels = c("Sinónimo (S)", "No sinónimo (N)", "Indel", "Codón de parada distinto", "Promotor", "UTR", "CDS sin cambios", "Intrón"),
                      name = NULL) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1B2A33"), guide = "none") +
    scale_x_continuous(limits = c(-9, xmax + 0.5), expand = c(0, 0)) +
    facet_wrap(~tramo, ncol = 1) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 8c. Mapa de todos los cambios entre especies en ", datos$gen),
         subtitle = paste0("Solo se muestran los codones que cambian respecto a ", nombre_esp(ref),
                           " (números = codón de ", nombre_esp(ref), "); el resto del gen está resumido con el mismo código de color.",
                           if (omitidos) paste0("\n", omitidos, " codones con cambios solo sinónimos no se muestran (aumenta max_codones para verlos).") else "")) +
    tema_base(11) + theme(axis.text = element_blank(), legend.position = "bottom",
                          strip.text = element_text(hjust = 0, face = "bold"))
  mostrar(p)
  guardar_figura(p, "fig8c_mapa_cambios", min(24, 3 + xmax * 0.13), 2.4 + ntr * (0.55 * ne * 2 + 0.7))
  resumen <- do.call(rbind, lapply(setdiff(esp, ref), function(e)
    data.frame(Especie = nombre_esp(e), S = sum(tipo[e, ] == "S"), N = sum(tipo[e, ] == "N"),
               Indel = sum(tipo[e, ] == "Indel"), `N/S` = round(sum(tipo[e, ] == "N") / max(1, sum(tipo[e, ] == "S")), 2),
               check.names = FALSE)))
  cat("  Codones con cambios respecto a", nombre_esp(ref), "(toda la CDS):\n")
  print(resumen, row.names = FALSE); cat("\n")
  datos$resumen_cambios <- resumen
  preguntas("¿En qué exón se concentran los cambios no sinónimos (N)? ¿Y los sinónimos (S)?",
            "¿La razón N/S aumenta o disminuye con la distancia evolutiva? ¿Qué indica una razón N/S baja?",
            "¿Hay codones que cambian en todas las especies? ¿Qué podría significar?")
  invisible(datos)
}

# ---------- Notación HGVS (c.) -> mutación -----------------------------------
parsear_hgvs <- function(h) {
  h0 <- h; h <- sub("^c\\.", "", gsub("\\s", "", h))
  m <- regmatches(h, regexec("^([0-9]+)([ACGT])>([ACGT])$", h))[[1]]
  if (length(m)) return(list(tipo = "sustitucion", posicion = as.integer(m[2]), ref = m[3], nuevo = m[4], hgvs = h0))
  m <- regmatches(h, regexec("^([0-9]+)(_([0-9]+))?delins([ACGT]+)$", h))[[1]]
  if (length(m)) { a <- as.integer(m[2]); b <- if (nzchar(m[4])) as.integer(m[4]) else a
    return(list(tipo = "delins", posicion = a, longitud = b - a + 1, nuevo = m[5], hgvs = h0)) }
  m <- regmatches(h, regexec("^([0-9]+)(_([0-9]+))?del([ACGT]*)$", h))[[1]]
  if (length(m)) { a <- as.integer(m[2]); b <- if (nzchar(m[4])) as.integer(m[4]) else a
    return(list(tipo = "delecion", posicion = a, longitud = b - a + 1, ref = if (nzchar(m[5])) m[5] else NULL, hgvs = h0)) }
  m <- regmatches(h, regexec("^([0-9]+)(_([0-9]+))?dup([ACGT]*)$", h))[[1]]
  if (length(m)) { a <- as.integer(m[2]); b <- if (nzchar(m[4])) as.integer(m[4]) else a
    return(list(tipo = "duplicacion", posicion = a, longitud = b - a + 1, ref = if (nzchar(m[5])) m[5] else NULL, hgvs = h0)) }
  m <- regmatches(h, regexec("^([0-9]+)_([0-9]+)ins([ACGT]+)$", h))[[1]]
  if (length(m)) return(list(tipo = "insercion", posicion = as.integer(m[2]), nuevo = m[4], hgvs = h0))
  stop("No reconozco la notación '", h0, "'. Ejemplos válidos: c.20A>T, c.68_69del, c.5266dup, c.30_31insGGG, c.10_12delinsAA")
}

# ---------- MÓDULO 7: estructura 3D de la proteína ---------------------------
# Usa el modelo de AlphaFold (EMBL-EBI) de la proteína y lo muestra con r3dmol:
# se puede ROTAR con el mouse (panel Viewer de RStudio o en el navegador).
# También guarda el archivo .pdb y scripts para PyMOL y ChimeraX.

COL_PLDDT <- c("Muy alta (>90)" = "#0053D6", "Alta (70-90)" = "#65CBF3",
               "Baja (50-70)" = "#FFDB13", "Muy baja (<50)" = "#FF7D45")

buscar_uniprot <- function(simbolo, taxon, proteina) {
  u <- sprintf("https://rest.uniprot.org/uniprotkb/search?query=%%28gene_exact%%3A%s%%29%%20AND%%20%%28organism_id%%3A%s%%29&fields=accession%%2Creviewed%%2Clength%%2Csequence&format=tsv&size=10",
               utils::URLencode(simbolo, reserved = TRUE), taxon)
  txt <- tryCatch(readLines(u, warn = FALSE), error = function(e) character())
  if (length(txt) < 2) return(NULL)
  t <- utils::read.delim(text = paste(txt, collapse = "\n"), stringsAsFactors = FALSE)
  if (!nrow(t)) return(NULL)
  k <- which(t$Sequence == proteina)[1]
  if (is.na(k)) k <- which(t$Reviewed == "reviewed")[1]
  if (is.na(k)) k <- which.min(abs(t$Length - nchar(proteina)))
  list(acc = t$Entry[k], identica = identical(t$Sequence[k], proteina), revisada = t$Reviewed[k] == "reviewed")
}

descargar_alphafold <- function(acc, destino) {
  info <- tryCatch(jsonlite::fromJSON(paste0("https://alphafold.ebi.ac.uk/api/prediction/", acc)),
                   error = function(e) NULL)
  url_pdb <- if (!is.null(info) && length(info$pdbUrl)) info$pdbUrl[1] else
    paste0("https://alphafold.ebi.ac.uk/files/AF-", acc, "-F1-model_v4.pdb")
  ok <- tryCatch({ utils::download.file(url_pdb, destino, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  ok && file.exists(destino) && file.size(destino) > 1000
}

leer_pdb_ca <- function(archivo) {
  l <- readLines(archivo, warn = FALSE)
  l <- l[startsWith(l, "ATOM") & substr(l, 13, 16) == " CA "]
  cad <- substr(l, 22, 22)
  l <- l[cad == cad[1]]                       # primera cadena
  tres <- trimws(substr(l, 18, 20))
  uno <- names(Biostrings::AMINO_ACID_CODE)[match(tools::toTitleCase(tolower(tres)), Biostrings::AMINO_ACID_CODE)]
  data.frame(resi = as.integer(substr(l, 23, 26)), aa = ifelse(is.na(uno), "X", uno),
             b = as.numeric(substr(l, 61, 66)), cadena = cad[1], stringsAsFactors = FALSE)
}

# Relaciona cada residuo de nuestra proteína con un residuo del archivo PDB.
mapear_proteina_pdb <- function(proteina, ca) {
  sec_pdb <- paste(ca$aa, collapse = "")
  if (sec_pdb == proteina) return(ca$resi)
  aln <- alinear(AAStringSet(c(prot = proteina, pdb = sec_pdb)))
  a <- letras(aln[["prot"]]); b <- letras(aln[["pdb"]])
  ip <- cumsum(a != "-"); ib <- cumsum(b != "-")
  mapa <- rep(NA_integer_, nchar(proteina))
  sel <- a != "-" & b != "-"
  mapa[ip[sel]] <- ca$resi[ib[sel]]
  mapa
}

modulo7_estructura3d <- function(datos, especie = NULL, comparar_con = NULL,
                                 color = c("propiedades", "confianza", "diferencias"),
                                 uniprot = NULL, archivo_pdb = NULL, girar = FALSE) {
  color <- match.arg(color)
  especie <- especie %||% datos$ref
  m <- datos$modelos[[especie]]
  if (is.null(m)) stop("La especie '", especie, "' no está entre las cargadas.")
  caja(paste0("MÓDULO 7 · Estructura 3D de la proteína ", datos$gen, " (", nombre_esp(especie), ")"),
       "Se usa el modelo predicho por AlphaFold (base de datos del EMBL-EBI).",
       "Los aminoácidos se colorean por propiedades químicas y se marcan los sitios",
       "de la mutación (módulo 6) y las diferencias con otra especie.",
       "Gira la proteína con el mouse; con la rueda acercas o alejas.")
  for (pq in c("r3dmol", "jsonlite", "htmlwidgets"))
    if (!requireNamespace(pq, quietly = TRUE)) install.packages(pq)
  dir.create(ruta("estructuras"), showWarnings = FALSE)

  # 1. Archivo de estructura
  if (is.null(archivo_pdb)) {
    info <- ESPECIES_INFO[ESPECIES_INFO$clave == especie, ]
    acc <- uniprot
    if (is.null(acc)) {
      cat("  Buscando la proteína en UniProt... ")
      u <- buscar_uniprot(m$simbolo, info$taxon, m$proteina)
      if (is.null(u)) stop("No se encontró la proteína en UniProt. Puedes indicar el código con uniprot = 'P68871' ",
                           "o un archivo propio con archivo_pdb = 'ruta/al/archivo.pdb'.")
      acc <- u$acc
      cat(acc, if (u$identica) "(secuencia idéntica a la del NCBI)" else "(secuencia similar; se alineará)", "\n")
    }
    archivo_pdb <- ruta("estructuras", paste0(datos$gen, "_", especie, "_AF-", acc, ".pdb"))
    if (!file.exists(archivo_pdb)) {
      cat("  Descargando el modelo de AlphaFold... ")
      if (!descargar_alphafold(acc, archivo_pdb)) stop("No se pudo descargar el modelo de AlphaFold para ", acc,
                                                     ". Revisa en ", url_alphafold(m$simbolo))
      cat("listo\n")
    }
  }
  ca <- leer_pdb_ca(archivo_pdb)
  mapa <- mapear_proteina_pdb(m$proteina, ca)
  cat(sprintf("  Estructura: %s (%d residuos; %d de %d aminoácidos de la proteína ubicados)\n",
              basename(archivo_pdb), nrow(ca), sum(!is.na(mapa)), nchar(m$proteina)))

  # 2. Clasificación de residuos
  aa <- letras(m$proteina)
  grupo <- GRUPOS_AA[aa]
  marcas <- list()
  if (!is.null(datos$mutacion) && datos$mutacion$especie == especie && !is.na(mapa[datos$mutacion$codon]))
    marcas$mutacion <- list(resi = mapa[datos$mutacion$codon], etiqueta = datos$mutacion$etiqueta_aa)
  difs <- integer()
  if (!is.null(comparar_con)) {
    if (is.null(datos$aln)) stop("Para comparar especies ejecuta antes el MÓDULO 1.")
    P <- matriz_aln(datos$aln$prot[c(especie, comparar_con)])
    pe <- columnas_a_posicion(P[1, ])
    dcol <- which(P[1, ] != P[2, ] & P[1, ] != "-")
    difs <- pe[dcol]
    cat(sprintf("  %d aminoácidos de %s difieren en %s.\n", length(difs), nombre_esp(especie), nombre_esp(comparar_con)))
  }

  # 3. Visor 3D rotable
  w <- r3dmol::r3dmol(backgroundColor = "white") |>
    r3dmol::m_add_model(data = archivo_pdb, format = "pdb") |>
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "#D9D9D9"))
  if (color == "propiedades") {
    for (gname in names(COL_GRUPOS)) {
      r <- mapa[which(grupo == gname)]; r <- r[!is.na(r)]
      if (length(r)) w <- w |> r3dmol::m_set_style(sel = r3dmol::m_sel(resi = r), style = r3dmol::m_style_cartoon(color = COL_GRUPOS[[gname]]))
    }
  } else if (color == "confianza") {
    corte <- cut(ca$b, c(-Inf, 50, 70, 90, Inf), labels = rev(names(COL_PLDDT)))
    for (nv in names(COL_PLDDT)) {
      r <- ca$resi[corte == nv]
      if (length(r)) w <- w |> r3dmol::m_set_style(sel = r3dmol::m_sel(resi = r), style = r3dmol::m_style_cartoon(color = COL_PLDDT[[nv]]))
    }
  } else if (color == "diferencias" && length(difs)) {
    r <- mapa[difs]; r <- r[!is.na(r)]
    w <- w |> r3dmol::m_set_style(sel = r3dmol::m_sel(resi = r), style = r3dmol::m_style_cartoon(color = "#E8833A"))
  }
  if (length(difs)) {
    r <- mapa[difs]; r <- r[!is.na(r)]
    if (length(r)) w <- w |> r3dmol::m_add_style(sel = r3dmol::m_sel(resi = r), style = r3dmol::m_style_stick(color = "#1B2A33", radius = 0.25))
  }
  if (!is.null(marcas$mutacion)) {
    w <- w |>
      r3dmol::m_add_style(sel = r3dmol::m_sel(resi = marcas$mutacion$resi), style = r3dmol::m_style_sphere(color = "#D62728")) |>
      r3dmol::m_add_label(text = marcas$mutacion$etiqueta, sel = r3dmol::m_sel(resi = marcas$mutacion$resi, atom = "CA"),
                          style = r3dmol::m_style_label(backgroundColor = "#D62728", fontSize = 14))
  }
  w <- w |> r3dmol::m_zoom_to()
  if (girar) w <- w |> r3dmol::m_spin()
  html <- ruta("estructuras", paste0(datos$gen, "_", especie, "_3D_", color, ".html"))
  autocontenido <- requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
  try(htmlwidgets::saveWidget(w, html, selfcontained = autocontenido, title = paste(datos$gen, "3D")), silent = TRUE)
  print(w)
  cat("  \U0001F9EC Visor 3D guardado: estructuras/", basename(html), " (ábrelo en el navegador y gíralo con el mouse)\n", sep = "")

  # 4. Scripts para PyMOL y ChimeraX
  hexrgb <- function(h) paste0("[", paste(round(grDevices::col2rgb(h)[, 1] / 255, 3), collapse = ", "), "]")
  pml <- c(paste0("load ", basename(archivo_pdb), ", prot"), "hide everything", "show cartoon", "color grey80")
  cxc <- c(paste0("open ", basename(archivo_pdb)), "hide atoms", "cartoon", "color #1 lightgray")
  for (gname in names(COL_GRUPOS)) {
    r <- mapa[which(grupo == gname)]; r <- r[!is.na(r)]
    if (!length(r)) next
    nom <- paste0("g_", gsub("[^A-Za-z]", "", iconv(gname, to = "ASCII//TRANSLIT")))
    pml <- c(pml, sprintf("set_color %s, %s", nom, hexrgb(COL_GRUPOS[[gname]])),
             sprintf("color %s, resi %s", nom, paste(r, collapse = "+")))
    cxc <- c(cxc, sprintf("color #1:%s %s", paste(r, collapse = ","), COL_GRUPOS[[gname]]))
  }
  if (length(difs)) {
    r <- mapa[difs]; r <- r[!is.na(r)]
    if (length(r)) {
      pml <- c(pml, sprintf("show sticks, resi %s and not name N+C+O", paste(r, collapse = "+")))
      cxc <- c(cxc, sprintf("show #1:%s atoms", paste(r, collapse = ",")))
    }
  }
  if (!is.null(marcas$mutacion)) {
    r <- marcas$mutacion$resi
    pml <- c(pml, sprintf("show spheres, resi %d", r), sprintf("color red, resi %d", r),
             sprintf('label resi %d and name CA, "%s"', r, marcas$mutacion$etiqueta))
    cxc <- c(cxc, sprintf("show #1:%d atoms", r), sprintf("style #1:%d sphere", r), sprintf("color #1:%d red", r),
             sprintf('label #1:%d text "%s"', r, marcas$mutacion$etiqueta))
  }
  pml <- c(pml, "orient")
  base <- ruta("estructuras", paste0(datos$gen, "_", especie))
  writeLines(pml, paste0(base, "_pymol.pml")); writeLines(cxc, paste0(base, "_chimerax.cxc"))
  cat("  Archivos para otros programas (misma carpeta 'estructuras'):\n",
      "   -", basename(archivo_pdb), "-> abre en Mol* (https://molstar.org/viewer), PyMOL o ChimeraX\n",
      "   -", basename(paste0(base, "_pymol.pml")), "-> en PyMOL: File > Run Script\n",
      "   -", basename(paste0(base, "_chimerax.cxc")), "-> en ChimeraX: open archivo.cxc\n\n")

  # 5. Figura 11: mapa lineal de la proteína (respaldo estático del 3D)
  n <- length(aa); ancho <- 50
  d <- data.frame(i = seq_len(n), aa = aa, grupo = factor(grupo, levels = names(COL_GRUPOS)))
  d$fila <- (d$i - 1) %/% ancho; d$col <- (d$i - 1) %% ancho + 1
  d$plddt <- ca$b[match(mapa, ca$resi)]
  d$conf <- cut(d$plddt, c(-Inf, 50, 70, 90, Inf), labels = rev(names(COL_PLDDT)))
  d$dif <- d$i %in% difs
  d$mut <- !is.null(datos$mutacion) && datos$mutacion$especie == especie
  d$mut <- d$mut & d$i == (datos$mutacion$codon %||% -1)
  H <- 2.2
  p11 <- ggplot(d) +
    geom_tile(aes(col, -fila * H, fill = grupo), colour = "white", height = 0.9) +
    geom_text(aes(col, -fila * H, label = aa), size = 2.8) +
    geom_tile(aes(col, -fila * H - 0.75, fill = conf), height = 0.35, colour = NA) +
    geom_point(data = d[d$dif, ], aes(col, -fila * H + 0.62), shape = 25, fill = "#1B2A33", size = 1.8) +
    geom_tile(data = d[d$mut, ], aes(col, -fila * H), fill = NA, colour = "#D62728", linewidth = 1.3, height = 0.95) +
    geom_text(data = data.frame(f = unique(d$fila)), aes(x = -1, y = -f * H, label = f * ancho + 1), hjust = 1, size = 3, colour = "#555555") +
    scale_fill_manual(values = c(COL_GRUPOS, COL_PLDDT), breaks = c(names(COL_GRUPOS), names(COL_PLDDT)), na.value = "#D9D9D9",
                      name = "Arriba: tipo de aminoácido · Abajo: confianza de AlphaFold (pLDDT)") +
    scale_x_continuous(limits = c(-4, ancho + 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         title = paste0("Figura 11. Proteína ", datos$gen, " (", nombre_esp(especie), "): propiedades, confianza de la estructura y sitios marcados"),
         subtitle = paste0("Recuadro rojo = sitio de la mutación", if (length(difs)) paste0(" · ▼ = aminoácido distinto en ", nombre_esp(comparar_con)) else "", ".")) +
    tema_base(11) + theme(axis.text = element_blank(), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
  mostrar(p11)
  guardar_figura(p11, paste0("fig11_mapa_proteina_", especie), 14, 2 + 0.55 * (max(d$fila) + 1))
  conf_media <- mean(d$plddt, na.rm = TRUE)
  if (grepl("AF-", basename(archivo_pdb))) cat(sprintf("  Confianza promedio del modelo (pLDDT): %.1f / 100\n\n", conf_media))
  else cat("  Nota: el archivo no es de AlphaFold; la franja inferior muestra el factor B, no la confianza pLDDT.\n\n")
  preguntas(
    "¿Los aminoácidos hidrofóbicos (naranja) quedan hacia el interior o la superficie de la proteína?",
    "¿El sitio de la mutación está en una hélice, una lámina o un bucle? ¿Está expuesto o escondido?",
    "¿Las diferencias entre especies (▼ y palitos negros) se concentran en zonas de baja confianza o en la superficie?",
    "¿Qué regiones tienen baja confianza (pLDDT < 50)? ¿Podrían ser regiones desordenadas?")
  preguntas_busqueda(
    paste0("Busca la proteína en UniProt (", url_uniprot(m$simbolo), "): en 'Family & Domains', ¿qué dominios tiene? ",
           "¿La mutación o las diferencias caen dentro de un dominio?"),
    paste0("En AlphaFold DB (", url_alphafold(m$simbolo), ") compara el modelo con una estructura experimental (PDB) si existe."))
  invisible(list(widget = w, mapa = mapa, archivo = archivo_pdb))
}
