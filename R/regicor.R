#' Càlcul de la puntuació regicor
#'
#' Funció vectoritzada per calcular els valors de l'escala de risc regicor
#'
#' @param df un dataframe amb els valors
#' @param noms_cols vector amb els noms de les columnes que s'utilitzaran pel càlcul. Vegeu detalls.
#'
#' @return Retorna un vector amb els valors de risc de regicor
#'
#' @keywords risc cardiovascular, regicor
#'
#' @details El paràmetre noms_cols és un vector amb els noms de les columnes del teu dataframe que s'utilitzaran
#' per fer els càlculs. S'han d'afegir els noms de les columnes per aquest ordre: edat, sexe, fumar, colesterol total,
#' tensió sistòlica, tensió diastòlica, colesterol HDL i diabetis.
#'
#' @export calcular_regicor
#'


calcular_regicor <- function(df, noms_cols) {
  noms_dict <- list(
    noms_cols[1] <- "Age",
    noms_cols[2] <- "Sex",
    noms_cols[3] <- "Smoking",
    noms_cols[4] <- "Total_cholesterol",
    noms_cols[5] <- "sBP",
    noms_cols[6] <- "dBP",
    noms_cols[7] <- "HDL_cholesterol",
    noms_cols[8] <- "DM"
  )

  for (i in 1:ncol(df)) {
    quin <- which(names(df)[i] %in% noms_cols)
    if (length(quin)) {
      names(df[, i]) <- noms_dict[quin]
    }
  }

  # valors dels punts de tall
  talls_col <- c(160, 200, 240, 280, Inf)
  talls_pas <- c(120, 130, 140, 160, Inf)
  talls_pad <- c(80, 85, 90, 100, Inf)
  talls_hdl <- c(35, 45, 50, 60, Inf)

  # valors dels coeficients
  # coeficients dones
  coefs_col_dones <- c(-0.26138, 0, 0.20771, 0.24385, 0.53513)
  coefs_ta_dones <- c(-0.53363, -0.06773, 0, 0.26288, 0.46573)
  coefs_hdl_dones <- c(0.84312, 0.37796, 0.19785, 0., -0.42951)

  # coeficients homes
  coefs_col_homes <- c(-0.65945, 0, 0.17692, 0.50539, 0.65713)
  coefs_ta_homes <- c(-0.00226, 0, 0.2832, 0.52168, 0.61859)
  coefs_hdl_homes <- c(0.49744, 0.2431, 0, -0.05107, -0.4866)

  # funcions

  # coeficients pel colesterol
  coefs_col <- function(df, talls_col, coefs_col_homes, coefs_col_dones) {
    temp <- matrix(unlist(lapply(df$Total_cholesterol, function(x) x < talls_col)), ncol = 5, byrow = T)
    vals <- apply(temp, 1, function(x) min(which(x)))
    ifelse(df$Sex == 1, coefs_col_homes[vals], coefs_col_dones[vals])
  }

  # coeficients per la TA
  coefs_ta <- function(df, talls_pas, talls_pad, coefs_ta_homes, coefs_ta_dones) {
    temp <- matrix(unlist(lapply(df$sBP, function(x) x < talls_pas)), ncol = 5, byrow = T)
    vals_sis <- apply(temp, 1, function(x) min(which(x)))
    temp <- matrix(unlist(lapply(df$dBP, function(x) x < talls_pad)), ncol = 5, byrow = T)
    vals_dia <- apply(temp, 1, function(x) min(which(x)))

    # Estem agafant el màxim de tots dos
    vals <- pmax.int(vals_dia, vals_sis)
    ifelse(df$Sex == 1, coefs_ta_homes[vals], coefs_ta_dones[vals])
  }

  # coeficients per hdl
  coefs_hdl <- function(df, talls_hdl, coefs_hdl_homes, coefs_hdl_dones) {
    temp <- matrix(unlist(lapply(df$HDL_cholesterol, function(x) x < talls_hdl)), ncol = 5, byrow = T)
    vals_hdl <- apply(temp, 1, function(x) min(which(x)))
    ifelse(df$Sex == 1, coefs_hdl_homes[vals_hdl], coefs_hdl_dones[vals_hdl])
  }

  # edat
  coef_edat <- function(df) {
    df$Age * ifelse(df$Sex == 1, 0.04826, 0.33766) + df$Age^2 * ifelse(df$Sex == 1, 0, -0.00268)
  }

  # diabets
  coefs_dm <- function(df) {
    ifelse(df$DM == 1, ifelse(df$Sex == 1, 0.42839, 0.59626), 0)
  }

  # fumar
  coefs_fum <- function(df) {
    ifelse(df$Smoking == 1, ifelse(df$Sex == 1, 0.52337, 0.29246), 0)
  }

  # càlcul del risc

  # càlcul del pre-risc
  pre_risc <- coef_edat(df) +
    coefs_col(df, talls_col, coefs_col_homes, coefs_col_dones) +
    coefs_hdl(df, talls_hdl, coefs_hdl_homes, coefs_hdl_dones) +
    coefs_ta(df, talls_pas, talls_pad, coefs_ta_homes, coefs_ta_dones) +
    coefs_dm(df) +
    coefs_fum(df)

  # càlcul de la probabilitat d'esdeveniment
  ifelse(df$Sex == 1, 1 - 0.951^exp(pre_risc - 3.489), 1 - 0.978^exp(pre_risc - 10.279))
}
