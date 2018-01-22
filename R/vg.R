####Identify high variability genes using Brennecke's method#####
findVG<-function(dat, plot.fig=NULL) {
        require(Matrix)
        require(ggplot2)
        if(!is.matrix(dat)){
          sf =Matrix::colSums(dat)
        }
        else{
          sf = colSums(dat)
        }
        sf = sf/median(sf)
	tmp.dat <- t( t(dat) / sf)
        if(is.matrix(dat)){
          g.means <- rowMeans(tmp.dat)
          g.vars <-  rowMeans(tmp.dat^2) - g.means^2
        }
        else{
          g.means <- Matrix::rowMeans(tmp.dat)
          g.vars <- Matrix::rowMeans(tmp.dat^2) - g.means^2
        }
	dispersion <- log10(g.vars/g.means)
        g.df = data.frame(g.means, g.vars, dispersion)
	#####test samples####
        dispersion = dispersion[!is.na(dispersion)]
        ###fit normal with 25% to 75%
        IQR = quantile(dispersion, c(0.25, 0.75))
        m = mean(IQR)
        delta = (IQR[2]-IQR[1])/(qnorm(0.75)-qnorm(0.25))
        z = (dispersion  - m)/delta        
        vg = data.frame(z)
        vg$pval = 1 - pnorm(z)
        vg$padj = p.adjust( vg$pval,method="fdr")
        vg = cbind(vg, g.df[row.names(vg),])

        ###loess regression
        select=vg$dispersion > 0
        fit=with(vg, loess(dispersion~log10(g.means), subset=select))
        rs=resid(fit)
        base=min(predict(fit))
        diff = vg$dispersion - base
        diff[select]  = rs
        IQR = quantile(diff, c(0.25, 0.75))
        m = mean(IQR)
        delta = (IQR[2]-IQR[1])/(qnorm(0.75)-qnorm(0.25))
        vg$loess.z = (diff  - m)/delta
        vg$loess.pval = 1 - pnorm(vg$loess.z)
        vg$loess.padj = p.adjust(vg$loess.pval, method="fdr")
        if(!is.null(plot.fig)){
          pdf(plot.fig)
          
          qqnorm(z)
          tmp=qqline(z)
          qqnorm(vg$loess.z)
          tmp=qqline(vg$loess.z)
          p=ggplot(vg, aes(z))+geom_density()
          plot(p)
          p=ggplot(vg, aes(loess.z))+geom_density()
          fit.df = data.frame(x=fit$x[,1],y=fit$fitted)
          plot(p)
          fit.df = fit.df[order(fit.df$x),]
          p=ggplot(vg, aes(x=log10(g.means), y=dispersion)) + geom_point() + geom_line(data=fit.df, aes(x=x,y=y,color="blue"))
          plot(p)
          dev.off()
        }
        return(vg)
      }

