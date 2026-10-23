# PsychoMatic 1.0.0: pruebas reproducibles del flujo simplificado.
# Ejecutar cada seccion por separado en RStudio. No requiere print()/summary().
library(PsychoMatic)
datos <- lavaan::HolzingerSwineford1939
items <- datos[paste0("x", 1:9)]
modelo <- 'visual =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed =~ x7 + x8 + x9'

desc_auto(items, digits = 4)
cormat(items, type = "pearson", digits = 4)
efa_auto(items, seed = 916, refine = "suggest")
cfa_auto(datos, modelo)  # seleccion automatica del estimador
cfa_auto(datos, modelo, estimator = "MLR", fit_variant = "scaled")
cfa_auto(datos, modelo, estimator = "MLR", fit_variant = "robust")
reliability_auto(items[1:3])
factorial_invariance_auto(datos, "school", modelo)

# Confiabilidad basada en el CFA: el objeto solo se guarda para reutilizarlo.
# El CFA y la confiabilidad imprimen sus resultados al ejecutarse.
ajuste <- cfa_auto(datos, modelo, estimator = "MLR")
reliability_auto(items, model = ajuste$fit)

# Exportar desde la propia funcion; los originales existentes se protegen.
carpeta <- tempfile("PsychoMatic-pruebas-")
dir.create(carpeta)
efa_auto(items, seed = 916, refine = "suggest", report = TRUE,
         format = c("excel", "word"), file_name = file.path(carpeta, "AFE"))
cfa_auto(datos, modelo, report = TRUE, format = "word",
         file_name = file.path(carpeta, "AFC.docx"))
reliability_auto(items[1:3], report = TRUE,
                 file_name = file.path(carpeta, "Confiabilidad.xlsx"))
factorial_invariance_auto(datos, "school", modelo, report = TRUE,
                         file_name = file.path(carpeta, "Invariancia.xlsx"))

# Opcionales: usar solo si responden a la pregunta de investigacion.
dimensionality_sensitivity(items, n_iter = 100, seed = 916)
dif_auto(items, group = datos$school)

# La invariancia parcial exige una liberacion y justificacion deliberadas.
# Ejemplo pedagogico, no recomendacion para modificar este modelo:
invariance_partial(datos, "school", modelo, free = "x2 ~ 1",
                   rationale = "Ejemplo pedagogico: evaluar sensibilidad del intercepto x2.",
                   estimator = "ML")

# Para procesos programaticos silenciosos:
resultado <- cfa_auto(datos, modelo, verbose = FALSE)

# Si una categoria ordinal esta vacia en un grupo, inspeccionar sus frecuencias:
# table(mis_datos$R5, mis_datos$grupo, useNA = "ifany")
# No eliminar ni fusionar categorias automaticamente para forzar un ajuste.
