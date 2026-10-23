# Validación de PsychoMatic

## Mantenimiento 1.0.1: corrección CRAN M1

Preparación local iniciada el 23 de octubre de 2026. CRAN ya publica 1.0.0;
su informe adicional M1mac identifica cuatro fallos de expectativas textuales
de rutas equivalentes `/var` y `/private/var`. Las correcciones canónicas que
ya estaban en main se incorporaron al nuevo tarball 1.0.1, con regresiones
adicionales de barras repetidas y creación/sobrescritura de reportes.

Los 44 motores estadísticos y las 47 fichas Rd siguen idénticos a CRAN 1.0.0.
La API no cambia. Los artefactos y PDF 1.0.0 se conservan como historial.
Archivo comprobado: `PsychoMatic_1.0.1.tar.gz`, 452715 bytes; SHA256
`30890BE9873CBBB2878FDA79B98C89D53AFF4C3D221D5339C32F08204D9535AD`.
El check completo `R CMD check --as-cran` en Windows 11 / R 4.6.1 finalizó
con **0 errores, 0 advertencias y 1 NOTE**: `Days since last update: 0`.
Esta NOTE administrativa se justifica porque CRAN solicitó corregir 1.0.0
el mismo día de su publicación; se conserva sin desactivar incoming.

Las pruebas aprobaron **3415 aserciones**, sin fallos ni advertencias, con
13 omisiones intensivas deliberadas del perfil de distribución. Ninguna de
las cuatro pruebas reportadas se omitió. Ejemplos, viñetas y manuales PDF y
HTML finalizaron OK. La revisión independiente no encontró bloqueos.

Manuales completos 1.0.1: EN 49 páginas, ES 51 páginas; 47 temas y 22 funciones.
SHA256 EN: `229C7E32AEC444DAB818E89A851B8547D1BF87B78C738F91C35C09AD3827EC57`.
SHA256 ES: `E44F2CD9B446CC6242D3970EA5DDF11715592E30AB44853624BAA62650DC9DCA`.

MacBuilder comprobó el mismo archivo en Apple M1 / macOS Tahoe 26.6 / R 4.6.1
Patched: **Status OK, 0 errores, 0 advertencias y 0 NOTEs**, 3416 aserciones
aprobadas y las mismas 13 omisiones previstas. La copia de los tests del
informe remoto coincide con la fuente comprobada. Consulte el
[informe de MacBuilder](https://mac.R-project.org/macbuilder/results/1792787197-53b23dc7692f1423/).
El servicio ejecutó su perfil R-release pese a solicitar R-devel mediante
el valor oficial del formulario; este resultado no acredita R-devel. Su check
estándar tampoco sustituye al `--as-cran` completo de Windows descrito arriba.
Los resultados del CI multiplataforma se consultan en
[GitHub Actions](https://github.com/gmoncayoj/PsychoMatic/actions/workflows/R-CMD-check.yaml)
y corresponden al commit fuente indicado en cada ejecución. Las comprobaciones
del tarball exacto se documentan por separado arriba. No se ha realizado el
reenvío editorial de 1.0.1 a CRAN.

## Historial: candidato GPL 1.0.0

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
