# Manuales completos de referencia / Complete reference manuals

Edición de mantenimiento 1.0.1: conserva las 22 funciones y 47 temas de 1.0.0;
la corrección afecta a las pruebas de rutas, sin cambios de API ni métodos.
Los manuales actualizados se construyeron y verificaron en ambos idiomas.

- [PsychoMatic_1.0.1_EN.pdf](PsychoMatic_1.0.1_EN.pdf): referencia inglesa actualizada.
- [PsychoMatic_1.0.1_ES.pdf](PsychoMatic_1.0.1_ES.pdf): traducción española actualizada.

Los PDF 1.0.0 se conservan para mantener válidos los enlaces del README que
ya forma parte del paquete publicado en CRAN.

Edición del 23 de octubre de 2026: ambos manuales incluyen GPL versión 3 o
posterior y copyright 2026 José Gamarra-Moncayo. El mantenedor conserva los PDFs
MIT previos como evidencia histórica local. Los hashes de los archivos vigentes
se encuentran en [VALIDATION.md](../VALIDATION.md).

- `PsychoMatic_1.0.0_EN.pdf`: referencia inglesa canónica.
- `PsychoMatic_1.0.0_ES.pdf`: traducción española completa.

Ambos incluyen todas las fichas de `man/`, también las que R marca como
internas: 47 temas, con las 22 funciones exportadas y métodos documentados.
Las firmas, opciones y sintaxis R son iguales en ambos idiomas. Los ejemplos
de ayuda usan plantillas que el investigador adapta a su data frame.

`man/` sigue siendo la ayuda inglesa instalada y la fuente que CRAN utiliza
para generar su manual oficial. El PDF español es documentación complementaria
del proyecto; no sustituye esa ayuda ni implica aprobación editorial de CRAN.
Los dos PDFs utilizan la composición nativa Rd/LaTeX del manual de referencia.

## Regenerar

Desde la raíz del repositorio, con R y una distribución LaTeX disponibles:

```sh
Rscript docs/manuals/check_manuals.R
Rscript docs/manuals/build_manuals.R
```

Actualice primero roxygen y la ayuda canónica cuando cambie la API; después
sincronice la prosa española de `es/man/`. `check_manuals.R` exige igualdad de
temas, nombres/alias, firmas, argumentos y ejemplos ejecutables, además de
cobertura de todos los exports. `coverage.csv` registra los hashes de ambas
fuentes. El constructor conserva sintaxis y etiquetas técnicas, traduce los
encabezados del manual español y usa una fecha explícita de edición.

Los archivos temporales se crean en `work/`, ignorado por Git. El directorio
`docs/` completo se excluye del paquete fuente por `.Rbuildignore`, evitando
incluir un segundo juego de ayuda o PDFs redundantes en el envío a CRAN.
