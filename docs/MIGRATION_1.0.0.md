# Migración de PsychoMatic 0.3.0 a 1.0.0

La serie de desarrollo 0.4.0 se integra en 1.0.0. Este número identifica la
primera API estable del proyecto y no significa que se haya publicado en CRAN.

## Código existente

Se conservan las funciones públicas previas y las posiciones de sus argumentos
compatibles. Revise los resultados antes de reutilizar un análisis archivado:
la selección de estimadores, preparación ordinal, diagnósticos y presentación
han recibido correcciones sustantivas documentadas en NEWS.

- En AFC e invarianza completa, `result$warnings` es una tabla estructurada.
  Para código antiguo que esperaba texto, utilice `result$warnings_list`, alias
  transitorio que se conserva durante el ciclo de deprecación.
- Los resultados nuevos suelen heredar de `psychomatic_result`: separan
  configuración solicitada/efectiva, resultados, diagnósticos y decisiones.
  La precisión de impresión no modifica los valores almacenados.
- La invarianza factorial y el procedimiento ordinal con identificación
  Wu–Estabrook tienen APIs diferentes. Seleccione el procedimiento de acuerdo
  con su modelo y explique cualquier liberación de parámetros en invarianza parcial.
- Las bases de demostración se mantienen como fixtures de desarrollo, fuera
  de la API instalada. Los ejemplos de ayuda son plantillas para los datos del
  investigador; las viñetas incluyen demostraciones reproducibles.
- Los informes automáticos se crean en `getwd()` con nombres únicos, salvo
  ruta explícita. `overwrite = TRUE` autoriza reemplazar archivos existentes.

## Funciones incorporadas

`psychomatic_spec()` declara estructura y codificación; `reliability_auto()`
centraliza fiabilidad; `ordinal_invariance_auto()` e `invariance_partial()`
amplían los procedimientos de invarianza; `network_auto()` y `network_plots()`
integran análisis de redes y figuras.

`esem_auto()` ajusta un bloque exploratorio de medición con lavaan. Para geomin,
seleccione columnas con `items = 2:13` o nombres, nunca pasando un data frame
en `items`. Target requiere una lista teórica de pertenencias primarias en
`factors`. Ambas rotaciones son oblicuas; el ajuste global no elige entre ellas.
`seed = NULL` conserva el comportamiento aleatorio del motor; una semilla
entera explícita permite reproducibilidad y restaura el estado previo del RNG.
La sección final de interpretación depende de los resultados y de su
admisibilidad; requiere revisión sustantiva antes de incorporarla a un manuscrito.

## Verificación antes de actualizar un proyecto

Guarde el código, la versión de R/dependencias y los resultados antiguos.
Instale primero el candidato en una biblioteca aislada, ejecute sus análisis
habituales y contraste preparación, muestras, estimadores e índices efectivos.
Los manuales de referencia completos EN/ES están en [manuals/](manuals/README.md)
y la evidencia resumida del candidato en [VALIDATION.md](VALIDATION.md).
