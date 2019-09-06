r::
install::
    sudo apt-get install r-base r-base-dev libopenblas-base

run::
    R

packages::
    install package in R, by
        install.packages("seismicRoll")

    load package in R, by
        library(seismicRoll)

data::
    load data in R, by
        data("iris");

function::
    define function, by
----------------------------------------------------------------------
        sta_lta <- function(x,ns,nl)
        {
            nx<-length(x) 
            sta<-array(NA,c(nx,1))
            lta<-array(NA,c(nx,1))
            
            for(i in 1:(nx-ns)){
                sta[i]=sum(x[i:(i+ns-1)])/ns
            }
            for(i in nl:nx){
                lta[i]=sum(x[(i-nl+1):i])/nl
            }
            r=sta/lta
            return(r)
        }
----------------------------------------------------------------------

    call function, by
        source("./sta_lta_r")
        x-<-sta_lta(x,ns,nl)
plot::
multiple lines::
    ylim should be set to let all graph has the same range

        plot(p,cex=1.5,col='red',type='b',ylim=c(0,max(x)))
        bline(v=first_break,col='red')
        par(new=T)
        plot(p1,type="l",ylim=c(0,max(x)))))"")
