# PsychoMatic

`PsychoMatic` es un paquete de R para automatizar analisis psicometricos y
facilitar su uso desde R y RStudio.

## Funciones incluidas

- `desc_auto()`: descriptivos y exportacion opcional a Excel.
- `efa_auto()`: analisis factorial exploratorio automatizado.
- `exportar_efa()`: exportacion de resultados de `efa_auto()` a Excel o Word.
- `cfa_auto()`: analisis factorial confirmatorio automatizado.
- `inv_align_auto()`: invarianza por alineamiento.
- `factorial_invariance_auto()`: invarianza factorial multigrupo.

## Instalacion desde GitHub

Sustituye `TU_USUARIO_GITHUB` por tu usuario real una vez subas el repositorio:

```r
install.packages("remotes")
remotes::install_github("TU_USUARIO_GITHUB/PsychoMatic")
```

## Flujo recomendado en RStudio

Abre la carpeta del paquete y ejecuta:

```r
devtools::document()
devtools::build()
devtools::check()
```

Si prefieres usar funciones base de R:

```r
roxygen2::roxygenise()
pkgbuild::build(".")
```

## Estructura del paquete

```text
PsychoMatic/
  DESCRIPTION
  NAMESPACE
  R/
  man/
  README.md
```

## Notas

- Actualiza el correo del mantenedor en `DESCRIPTION`.
- Cambia la licencia si quieres usar otra distinta.
- La documentacion y `NAMESPACE` ahora se regeneran con `roxygen2`.
