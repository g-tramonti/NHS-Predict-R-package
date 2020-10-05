#'@title scan.predict
#'
#'@description  Calculates 'NHS Predict' v2.1 Overall survival and chemotherapy benefits in a simplified version with fewer inputs, suited for SCAN audit.
#'
#'@param data A dataframe containing patient data with the necessary variables.Except for age at diagnosis, the other variables must be named according to SCAN
#             convention: MREFER, MAXPATH, TGRADE, INVNODES, ERSTATUS, HER2STATUS.
#'@param age.start Numeric, Age at diagnosis of the patient. Range between 25 and 85.
#'
#'@return The function attaches additional columns to the dataframe, matched for patient observation,
#'    containing Overall survival at the specified year, plus the additional benefit for chemotherapy at 5, 10, and 15 years interval.
#'    Observations containing missing values are moved to the bottom.
#'
#'@examples
#'data(scan_example_data)
#'
#'scan_example_data <- scan.predict(scan_example_data, age.start = diag_age)
#'
#'@importFrom stats complete.cases
#'
#'@export

scan.predict <- function(data,age.start){
  age.start  <- deparse(substitute(age.start))
  MREFER     <- deparse(substitute(MREFER))
  MAXPATH    <- deparse(substitute(MAXPATH))
  TGRADE     <- deparse(substitute(TGRADE))
  INVNODES   <- deparse(substitute(INVNODES))
  ERSTATUS   <- deparse(substitute(ERSTATUS))
  HER2STATUS <- deparse(substitute(HER2STATUS))



  predict2.1s <- function(age.start, MREFER, MAXPATH, TGRADE, INVNODES, ERSTATUS, HER2STATUS) {
    ## ---------------------------------------------------------------------
    #Fix inputs and default values
    age.start <- ifelse(is.na(age.start)|age.start<25|age.start>85, 1000, age.start)
    MREFER <- ifelse(MREFER==9,0.204,ifelse(MREFER!=2,0,1))
    TGRADE <- ifelse(TGRADE == 1 | TGRADE == 2 | TGRADE == 3, TGRADE,2.13)
    MAXPATH <- ifelse(MAXPATH>999|is.na(MAXPATH),1000,MAXPATH)
    INVNODES <- ifelse(INVNODES>999|is.na(INVNODES),1000,INVNODES)
    ERSTATUS <- ifelse(ERSTATUS == 1 | ERSTATUS == 3, 1, ifelse(ERSTATUS==2,0,99))
    HER2STATUS <- ifelse(HER2STATUS == 1, 1, ifelse(HER2STATUS == 2, 0, 9))
    ki67 <- 9
    generation <- 3
    horm <- ifelse(ERSTATUS == 1, 1, 0)
    ## ------------------------------------------------------------------------
    time      <- c(1:15)
    age       <- age.start - 1 + time
    grade.val <- ifelse(ERSTATUS==1, TGRADE, ifelse(TGRADE==2 || TGRADE == 3, 1, 0)) # Grade variable for ER neg
    ## ------------------------------------------------------------------------
    age.mfp.1   <- ifelse(ERSTATUS==1, (age.start/10)^-2-.0287449295, age.start-56.3254902)
    age.beta.1  <- ifelse(ERSTATUS==1, 34.53642, 0.0089827)
    age.mfp.2   <- ifelse(ERSTATUS==1, (age.start/10)^-2*log(age.start/10)-.0510121013, 0)
    age.beta.2  <- ifelse(ERSTATUS==1, -34.20342, 0)
    size.mfp    <- ifelse(ERSTATUS==1, log(MAXPATH/100)+1.545233938, (MAXPATH/100)^.5-.5090456276)
    size.beta   <- ifelse(ERSTATUS==1, 0.7530729, 2.093446)
    nodes.mfp   <- ifelse(ERSTATUS==1, log((INVNODES+1)/10)+1.387566896, log((INVNODES+1)/10)+1.086916249)
    nodes.beta  <- ifelse(ERSTATUS==1, 0.7060723, .6260541)
    grade.beta  <- ifelse(ERSTATUS==1, 0.746655, 1.129091)
    screen.beta <- ifelse(ERSTATUS==1, -0.22763366, 0)
    her2.beta   <- ifelse(HER2STATUS==1, 0.2413,
                          ifelse(HER2STATUS==0, -0.0762,0 ))
    ki67.beta   <- ifelse(ki67==1 & ERSTATUS==1, 0.14904,
                          ifelse(ki67==0 & ERSTATUS==1, -0.11333,0 ))

    ## ----baseline_adjust-----------------------------------------------------
    r.prop   <- 0.69 # Proportion of population receiving radiotherapy
    r.breast <- 0.82 # Relative hazard breast specific mortality
    r.other  <- 1.07 # Relative hazard other mortality

    ## ------------------------------------------------------------------------
    # Other mortality prognostic index (mi)
    mi <- 0.0698252*((age.start/10)^2-34.23391957)

    # Breast cancer mortality prognostic index (pi)
    pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
      nodes.beta*nodes.mfp + grade.beta*grade.val + screen.beta*MREFER +
      her2.beta + ki67.beta

    c     <- -.446
    h     <- ifelse(horm==1 & ERSTATUS==1, -0.3857, 0)

    hc <- h + c

    ## ------------------------------------------------------------------------
    # Generate cumulative baseline other mortality
    base.m.cum.oth <- exp(-6.052919 + (1.079863*log(time)) + (.3255321*time^.5))

    # Generate cumulative survival non-breast mortality
    s.cum.oth <- exp(-exp(mi)*base.m.cum.oth)

    # Generate annual baseline non-breast mortality
    base.m.oth <- base.m.cum.oth
    for (i in 2:15) {
      base.m.oth[i] <- base.m.cum.oth[i] - base.m.cum.oth[i-1] }

    # Loop for different treatment options
    rx.oth <- c(surg = 0, h = 0, c = 0, hc = 0)

    cols <- length(rx.oth) # Number of Treatment categories

    # Generate the annual non-breast mortality rate
    # Matrix (15x9) with column for each treatment
    m.oth.rx <- sapply(rx.oth, function(rx.oth, x.vector = base.m.oth) {
      output <-  x.vector*exp(mi + rx.oth)
      return(output)
    }
    )

    # Calculate the cumulative other mortality rate
    m.cum.oth.rx <- apply(m.oth.rx, 2, cumsum)

    # Calculate the cumulative other survival
    s.cum.oth.rx <- exp(-m.cum.oth.rx)

    # Convert cumulative mortality rate into cumulative risk
    m.cum.oth.rx <- 1- s.cum.oth.rx

    m.oth.rx <- m.cum.oth.rx
    for (j in 1:cols) {
      for (i in 2:15) {
        m.oth.rx[i,j] <- m.cum.oth.rx[i,j] - m.cum.oth.rx[i-1,j]
      }
    }

    ## ------------------------------------------------------------------------
    # Generate cumulative baseline breast mortality
    if (ERSTATUS==1) {
      base.m.cum.br <- exp(0.7424402 - 7.527762/time^.5 - 1.812513*log(time)/time^.5)
    } else { base.m.cum.br <- exp(-1.156036 + 0.4707332/time^2 - 3.51355/time)
    }

    # Generate annual baseline breast mortality
    base.m.br <- base.m.cum.br
    for (i in 2:15) {
      base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }

    # Loop for different treatment options
    rx <- c(surg = 0, h = h, c = c, hc = hc)

    # Generate the annual breast cancer specific mortality rate
    m.br.rx <- sapply(rx, function(rx, x.vector = base.m.br) {
      output <-  x.vector*exp(pi + rx)
      return(output)
    }
    )

    # Calculate the cumulative breast cancer mortality rate
    m.cum.br.rx <- apply(m.br.rx, 2, cumsum)

    # Calculate the cumulative breast cancer survival
    s.cum.br.rx <- exp(- m.cum.br.rx)

    # Convert cumulative mortality rate into cumulative risk
    m.cum.br.rx <- 1- s.cum.br.rx

    m.br.rx <- m.cum.br.rx
    for (j in 1:cols) {
      for (i in 2:15) {
        m.br.rx[i,j] <- m.cum.br.rx[i,j] - m.cum.br.rx[i-1,j]
      }
    }

    ## ------------------------------------------------------------------------
    m.cum.all.rx <- 1 - s.cum.oth.rx*s.cum.br.rx
    s.cum.all.rx <- 100-100*m.cum.all.rx

    # Annual all cause mortality
    m.all.rx <- m.cum.all.rx
    for (j in 1:cols) {
      for (i in 2:15) {
        m.all.rx[i,j] <- m.cum.all.rx[i,j] - m.cum.all.rx[i-1,j]
      }
    }

    ## ------------------------------------------------------------------------

    # Proportion of all cause mortality that is breast cancer
    prop.br.rx      <- m.br.rx/(m.br.rx + m.oth.rx)
    pred.m.br.rx    <- prop.br.rx*m.all.rx
    pred.cum.br.rx  <- apply(pred.m.br.rx, 2, cumsum)
    pred.m.oth.rx   <- m.all.rx - pred.m.br.rx
    pred.cum.oth.rx <- apply(pred.m.oth.rx, 2, cumsum)
    pred.cum.all.rx <- pred.cum.br.rx + pred.cum.oth.rx

    ## ------------------------------------------------------------------------
    # Treatment benefits

    benefits2.1 <- 100*(pred.cum.all.rx[,1] - pred.cum.all.rx)

    Overall     <- round((1-pred.cum.all.rx[10,1])*100 + benefits2.1[10,2],2) #Overall 10 years
    #chemo.ben5  <- round(benefits2.1[5,4]-benefits2.1[5,2],2) #Chemotherapy benefit 5 years
    chemo.ben10 <- round(benefits2.1[10,4]-benefits2.1[10,2],2) #Chemotherapy benefit 10 years
    #chemo.ben15 <- round(benefits2.1[15,4]-benefits2.1[15,2],2) #Chemotherapy benefit 15 years
    Overall <- ifelse(age.start==1000|INVNODES==1000|MAXPATH==1000|ERSTATUS==99,NA,Overall)
    chemo.ben10 <- ifelse(is.na(Overall),NA,chemo.ben10)
    return(c(Overall, chemo.ben10))
  }

  {N<-nrow(data) #Number of patients
    results <- array(c(NA,NA),c(N,2))
    O.S  <- paste("Overall.10years")
    B.10 <- paste("Chemo.Benefit.10years")
    colnames(results)<-c(O.S,B.10)} #Rename columns

  for (n in 1:N) {
    results[n,]<- predict2.1s(data[[age.start]][n],data[[MREFER]][n],data[[MAXPATH]][n],
                              data[[TGRADE]][n],data[[INVNODES]][n],data[[ERSTATUS]][n],data[[HER2STATUS]][n])
  }

  data <- cbind(data, results) #Paste results to data
  temp.dat <- complete.cases(data)
  data <- rbind(data[temp.dat,],data[!temp.dat,])
}
