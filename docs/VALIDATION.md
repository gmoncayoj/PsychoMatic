# Validación de PsychoMatic 1.0.0

Comprobaciones finalizadas el 23 de octubre de 2026. Este documento describe
el candidato con licencia GPL versión 3 o posterior y copyright 2026
José Gamarra-Moncayo. No constituye una confirmación de publicación en CRAN.

## Archivo comprobado

`PsychoMatic_1.0.0.tar.gz`, 459506 bytes.
SHA256: `A062B6F2C4107BA1D1207907A2DBE61C890E21DB259AB10C0C10E4C68BA5B9BC`.
El hash identifica el archivo comprobado; una reconstrucción puede producir
otro hash por sus metadatos, incluso conservando el código.

| Entorno | Resultado del mismo tarball GPL |
| --- | --- |
| Windows 11 x64, R 4.6.1 | `R CMD check --as-cran` completo, Status OK |
| Windows x64, R-devel r90748 (2026-10-22) | Win-builder, Status OK |

Ambas comprobaciones finalizaron sin errores, advertencias ni NOTEs e incluyeron
pruebas, ejemplos, viñetas y manuales PDF/HTML. El
[informe de Win-builder](https://win-builder.r-project.org/AyEWPJy5LMt4/00check.log)
documenta el resultado remoto; su copia y el informe local se conservan en la
evidencia del mantenedor porque los enlaces del servicio son temporales.

La suite instalada local aprobó **3402 aserciones**, con cero fallos y
advertencias y **13 omisiones intensivas deliberadas**. Esas regresiones siguen
disponibles con `PSYCHOMATIC_EXTENDED_TESTS=true`; consulte [TESTING.md](TESTING.md).
Las pruebas tardaron 583.53 segundos en esa ejecución. El tiempo varía con la
carga y el entorno; no es un benchmark ni una garantía de duración.

Los 44 archivos R, 47 fichas de ayuda y 60 archivos de pruebas del candidato
GPL son idénticos a los del candidato MIT precedente. El cambio de licencia se
comprobó dentro del nuevo archivo. Los checks Linux/macOS del candidato anterior
no se atribuyen al nuevo tarball GPL. La matriz de GitHub Actions está
configurada para Linux, Windows, macOS y distintas versiones de R, pero no se
había ejecutado en GitHub al cerrar estas comprobaciones.

## Manuales y dependencias inversas

Los [manuales completos](manuals/README.md) abarcan las 22 funciones públicas y
los 47 temas de ayuda; se verificaron las mismas firmas y ejemplos entre idiomas.

| Manual | Páginas | SHA256 |
| --- | --- | --- |
| Inglés | 49 | `53B5915329A98DF9464F123C1FC288CD4DA30AA078F393E8AFA6D6B566AE5148` |
| Español | 51 | `76B10102D01418C8369978D183FD35AF30413F36C5CBA85794FA7DAAEE197C5A` |

El índice CRAN revisado el 22 de octubre de 2026 contenía 25528 paquetes y no
mostró dependencias inversas de PsychoMatic en Depends, Imports, LinkingTo,
Suggests o Enhances, incluidas las recursivas. Esto no cubre código privado ni
repositorios externos. El tutorial de lavaan y sus fuentes oficiales también
respondieron HTTP200 con el contenido esperado en esa fecha.
