#'@title rfs.predict
#'
#'@description  Calculates Predict2.1 Recurrence-free survival and therapy benefits
#'
#'@param data A dataframe containg patient data with the necessary variables.
#'@param year Numeric, Specify the year since surgery for which the predictions are calculated, ranges between 1 and 15. Default at 10.
#'@param age.start Numeric, Age at diagnosis of the patient. Range between 25 and 85.
#'@param screen Numeric, Clinically detected = 0, Screen detected = 1.
#'@param size Numeric, Tumour size in millimiters.
#'@param grade Numeric, Tumour grade. Values: 1,2,3. Missing=9.
#'@param nodes Numeric, Number of positive nodes.
#'@param er Numeric, ER status, ER+ = 1, ER- = 0.
#'@param her2 Numeric, HER2 status, HER2+ = 1, HER2- = 0. Unknown = 9.
#'@param ki67 Numeric, ki67 status, KI67+ = 1, KI67- = 0, Unknown = 9.
#'@param generation Numeric, Chemotherapy generation. Values: 0,2,3.
#'@param horm Numeric, Hormone therapy, Yes = 1, No = 0.
#'@param traz Numeric, Trastuzumab therapy, Yes = 1, No = 0.
#'@param bis Numeric, Bisphosphonate therapy, Yes = 1, No = 0.
#'
#'@return The function attaches additional columns to the dataframe, matched for patient observation,
#'    containing recurrence-free survival at the specified year, plus the additional benefit for each type of therapy.
#'
#'@examples
#'
#'data(example_data)
#'
#'example_data <- rfs.predict(example_data,age.start = age,screen = detection,size = t.size,
#'                      grade = t.grade, nodes = nodes, er = er.status, her2 = her2.status,
#'                      ki67 = ki67.status, generation = chemo.gen, horm = horm.t,
#'                      traz = trastuzumab, bis = bis.t)
#'@examples
#'
#'data(example_data)
#'
#'example_data <- rfs.predict(example_data,year = 15, age,detection,t.size,t.grade,
#'                            nodes,er.status,her2.status,ki67.status,chemo.gen,horm.t,
#'                            trastuzumab,bis.t)
#'
#'@export

rfs.predict <- function(data,year=10,age.start,screen,size,grade,nodes,er,her2,ki67,generation,horm,
                        traz,bis){
  age.start <- deparse(substitute(age.start))
  screen <- deparse(substitute(screen))
  size <- deparse(substitute(size))
  grade <- deparse(substitute(grade))
  nodes <- deparse(substitute(nodes))
  er <- deparse(substitute(er))
  her2 <- deparse(substitute(her2))
  ki67 <- deparse(substitute(ki67))
  generation <- deparse(substitute(generation))
  horm <- deparse(substitute(horm))
  traz <- deparse(substitute(traz))
  bis <- deparse(substitute(bis))

  predict2.1.2 <- function(age.start,screen,size,grade,nodes,er,her2,ki67,generation,horm,
                          traz,bis){
    ##----------------------------------------------------------------
    #Fix inputs
    screen    <- ifelse(screen == 2, 0.204, screen)
    grade     <- ifelse(grade == 9, 2.13, grade)

    ## ------------------------------------------------------------------------
    time      <- c(1:15)
    age       <- age.start - 1 + time
    grade.val <- ifelse(er==1, grade, ifelse(grade==2 || grade == 3, 1, 0)) # Grade variable for ER neg
    ## ------------------------------------------------------------------------
    age.mfp.1   <- ifelse(er==1, (age.start/10)^-2-.0287449295, age.start-56.3254902)
    age.beta.1  <- ifelse(er==1, 34.53642, 0.0089827)
    age.mfp.2   <- ifelse(er==1, (age.start/10)^-2*log(age.start/10)-.0510121013, 0)
    age.beta.2  <- ifelse(er==1, -34.20342, 0)
    size.mfp    <- ifelse(er==1, log(size/100)+1.545233938, (size/100)^.5-.5090456276)
    size.beta   <- ifelse(er==1, 0.7530729, 2.093446)
    nodes.mfp   <- ifelse(er==1, log((nodes+1)/10)+1.387566896, log((nodes+1)/10)+1.086916249)
    nodes.beta  <- ifelse(er==1, 0.7060723, .6260541)
    grade.beta  <- ifelse(er==1, 0.746655, 1.129091)
    screen.beta <- ifelse(er==1, -0.22763366, 0)
    her2.beta   <- ifelse(her2==1, 0.2413,
                          ifelse(her2==0, -0.0762,0 ))
    ki67.beta   <- ifelse(ki67==1 & er==1, 0.14904,
                          ifelse(ki67==0 & er==1, -0.11333,0 ))

    ## ----baseline_adjust-----------------------------------------------------
    r.prop   <- 0.69 # Proportion of population receiving radiotherapy
    r.breast <- 0.82 # Relative hazard breast specific mortality
    r.other  <- 1.07 # Relative hazard other mortality

    ## ------------------------------------------------------------------------
    # Other mortality prognostic index (mi)
    mi <- 0.0698252*((age.start/10)^2-34.23391957)

    # Breast cancer mortality prognostic index (pi)
    pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
      nodes.beta*nodes.mfp + grade.beta*grade.val + screen.beta*screen +
      her2.beta + ki67.beta

    c     <- ifelse(generation == 0, 0, ifelse(generation == 2, -.248, -.446))
    h     <- ifelse(horm==1 & er==1, -0.3857, 0)
    t     <- ifelse(her2==1 & traz==1, -.3567, 0)
    b     <- ifelse(bis==1, -0.198, 0) # Only applicable to menopausal women.


    hc <- h + c
    ht <- h + t
    hb <- h + b
    ct <- c + t
    cb <- c + b
    tb <- t + b
    hct <- h + c + t
    hcb <- h + c + b
    htb <- h + t + b
    ctb <- c + t + b
    hctb <- h + c + t + b

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
    rx.oth <- c(surg = 0,
                h = 0,
                c = 0,
                t = 0,
                b = 0,
                hc = 0,
                ht = 0,
                hb = 0,
                ct = 0,
                cb = 0,
                tb = 0,
                hct = 0,
                hcb = 0,
                htb = 0,
                ctb = 0,
                hctb = 0)

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
    if (er==1) {
      base.m.cum.br <- exp(0.7424402 - 7.527762/time^.5 - 1.812513*log(time)/time^.5)
    } else { base.m.cum.br <- exp(-1.156036 + 0.4707332/time^2 - 3.51355/time)
    }

    # Loop for different treatment options
    rx <- c(surg = 0,
            h = h,
            c = c,
            t = t,
            b = b,
            hc = hc,
            ht = ht,
            hb = hb,
            ct = ct,
            cb = cb,
            tb = tb,
            hct = hct,
            hcb = hcb,
            htb = htb,
            ctb = ctb,
            hctb = hctb)

    # Generate annual baseline breast mortality
    base.m.br <- base.m.cum.br
    for (i in 2:15) {
      base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }

    dfs <- c(rep(0.687,5), rep(0.016,10))
    base.rel <- base.m.br*exp(dfs)

    # # Generate the baseline annual relapse rate
    base.rel.rx <- sapply(rx, function(rx, x.vector = base.rel) {
      output <-  x.vector*exp(pi + rx)
      return(output)
     }
    )

    # # Calculate the cumulative relapse rate
     rel.cum.rx <- apply(base.rel.rx, 2, cumsum)

    # # Calculate the cumulative relapse survival
     s.cum.rel.rx <- exp(- rel.cum.rx)

    # # Convert cumulative relapse rate into cumulative risk
     rel.cum.rx <- 1- s.cum.rel.rx

    # # Calculate annual relapse risk
     rel.rx <- rel.cum.rx
     for (j in 1:cols) {
       for (i in 2:15) {
         rel.rx[i,j] <- rel.cum.rx[i,j] - rel.cum.rx[i-1,j]
       }
     }

    # # Cumulative all cause survival conditional on surviving relapse and all cause mortality
     m.cum.all.rx <- 1 - s.cum.oth.rx*s.cum.rel.rx
     s.cum.all.rx <- 100-100*m.cum.all.rx

    # # Annual all cause event rate
     m.all.rx <- m.cum.all.rx
     for (j in 1:cols) {
       for (i in 2:15) {
         m.all.rx[i,j] <- m.cum.all.rx[i,j] - m.cum.all.rx[i-1,j]
       }
     }

    # # Proportion of all events that is relapse
     prop.rel.rx <- rel.rx/(rel.rx + m.oth.rx)

    # # Calulate the predicted relapse and other mortality
     pred.m.rel.rx       <- prop.rel.rx*m.all.rx
     pred.cum.rel.rx     <- apply(pred.m.rel.rx, 2, cumsum)
     pred.m.oth.rx       <- m.all.rx - pred.m.rel.rx
     pred.cum.oth.rx     <- apply(pred.m.oth.rx, 2, cumsum)
     pred.cum.all.rel.rx <- pred.cum.rel.rx + pred.cum.oth.rx

    # rx benefits

    benefits2.1.2 <- 100*(pred.cum.all.rel.rx[,1] - pred.cum.all.rel.rx)

    Overall <- round((1-pred.cum.all.rel.rx[year,1])*100,2) #Overall
    horm.ben <- round(benefits2.1.2[year,2],2) #Hormonal therapy benefit
    chemo.ben <- round(benefits2.1.2[year,6]-benefits2.1.2[year,2],2) #Chemotherapy benefit
    traz.ben <- round(benefits2.1.2[year,12]-benefits2.1.2[year,6],2) #Trastuzumab benefit
    bis.ben <- round(benefits2.1.2[year,16]-benefits2.1.2[year,12],2) #Bisphosphonate benefit

    return(c(Overall,horm.ben,chemo.ben,traz.ben,bis.ben))
  }

  {N<-nrow(data) #Number of patients
    results <- array(c(NA,NA),c(N,5))
    R.S <- paste("Relapse-free.",year,".years")
    H.B <- paste("Horm.Benefit.",year,".years")
    C.B <- paste("Chemo.Benefit.",year,".years")
    T.B <- paste("Traz.Benefit.",year,".years")
    B.B <- paste("Bis.Benefit.",year,".years")
    colnames(results)<-c(R.S,H.B,C.B,T.B,B.B)} #Rename columns

  for (n in 1:N) {
    results[n,]<- predict2.1.2(data[[age.start]][n],data[[screen]][n],data[[size]][n],data[[grade]][n],data[[nodes]][n],
                              data[[er]][n],data[[her2]][n],data[[ki67]][n],data[[generation]][n],data[[horm]][n],
                              data[[traz]][n],data[[bis]][n])
  }
  data <- cbind(data, results) #Paste results to data
}



