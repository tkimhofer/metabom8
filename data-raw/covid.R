## code to prepare `DATASET` dataset goes here

usethis::use_data("covid")

load('/Users/TKimhofer/Rproj/cov19/CoV_rawCPMG_calibrated.Rdata')
Xc[is.na(Xc)]=0
an=an[,colnames(an) %in% c('type', 'hospital')]
an$lineWidthHz=round(lw(Xc, ppm)*meta$a_SFO1,1)
an$noiseLevel=round(noise.est(Xc, ppm))

# calibrate to TSP (use examples to calibrate to glucose)
Xc=calibrate(Xc, ppm, type='glucose')
#matspec(Xc, ppm, shift=c(5.15, 5.3))

# cap ends and excise res. water signal
#matspec(Xc, ppm, shift=c(4.5, 4.9))

idx=c(get.idx(c(0.25, 4.55), ppm), get.idx(c(4.8, 9), ppm))
ppm=ppm[idx]
X=Xc[, idx]

idx=c(which(an$type=='Patients Cov19 (+)')[c(1:3, (length(which(an$type=='Patients Cov19 (+)'))-1):length(which(an$type=='Patients Cov19 (+)')))], sample(which(an$type=='Healthy'), 5))
an=an[idx,]
X=X[idx,]
meta=meta[idx,]

rownames(an)=1:nrow(an)
rownames(X)=rownames(an)
rownames(meta)=rownames(an)

matspec(X, ppm, shift=range(ppm))

# bl correction
Xbl=bline(X)

# normalise
Xn=pqn(Xbl, add_DilF = 'dilf')

X=Xn
save(X, ppm, an, file='data/covid.rda')
