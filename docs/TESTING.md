# Pruebas reproducibles

La suite de distribución conserva validación de entradas, decisiones, ajuste
nativo de modelos, equivalencia de resultados y reportes. Las regresiones
extensas se mantienen en el mismo árbol de pruebas y se habilitan explícitamente.
No se modifica ningún motor estadístico ni sus valores predeterminados para
acelerar las comprobaciones.

La selección extendida contiene las matrices exhaustivas de datos ausentes en
ESEM, una matriz extensa de variantes de invarianza parcial, soluciones nativas
impropias/no convergentes, exportaciones completas repetidas y las integraciones
de redes con más réplicas, repetición de semilla y todas las figuras. Sus cuerpos
se conservan. La selección ordinaria mantiene equivalencia nativa ESEM
target/geomin tanto continua como ordinal, invarianza, estados no evaluables,
reportes bilingües y el smoke nativo de redes. No se omiten fallos conocidos.

## Suite de distribución

```r
Sys.unsetenv("PSYCHOMATIC_EXTENDED_TESTS")
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```

El tarball ejecuta esta misma selección con `R CMD check --as-cran`. Las
omisiones identificadas como `Extended regression` son deliberadas y se
distinguen de omisiones por dependencias ausentes.

## Suite completa

```r
Sys.setenv(PSYCHOMATIC_EXTENDED_TESTS = "true")
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```

La opción debe ser exactamente `true`. No depende de `interactive()`, de
`NOT_CRAN`, de una versión de R o de la velocidad del equipo. El workflow
`.github/workflows/R-CMD-check.yaml` la establece explícitamente para conservar
las regresiones nativas completas en Linux, Windows y macOS. Configurar un
workflow no demuestra que se haya ejecutado en GitHub.

## Evidencia y límites

El [resumen de validación](VALIDATION.md) identifica el candidato GPL comprobado.
El mantenedor conserva localmente perfiles por caso, comprobación de la política
de selección, revisión de dependencias inversas, verificación HTTP y resultados
de los comprobadores. Cada plataforma debe identificar el SHA256 del tarball;
un envío pendiente no cuenta como aprobación. La revisión del índice CRAN no
detecta dependencias privadas o repositorios externos.

Los smoke compactos vuelven a abrir Word y Excel en ambos idiomas y comprueban
precisión numérica, y ejecutan bootstrap y permutaciones de redes con motores
nativos y pocas réplicas. Ese presupuesto pequeño sirve para verificar la
integración del software; no constituye una configuración recomendada para
analizar datos de investigación.
