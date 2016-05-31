
  prediction.surface.kriging = function( p, DS="annual", method="direct" , add.set.locations=F, vclass="R0.mass" ) {
    # prepare surface area and temperature estimates for plots
    
    if (DS %in% c("annual", "annual.redo" ) ) {

      fn = file.path( project.datadirectory("snowcrab"), "R", "PS.habitat.ts.rdata" )    
      
      if (DS=="annual") {
        load(fn)
        return( PS.hab )
      }
      
      out = NULL
      
      # yrs =  p$yearswithTdata --- no longer possible as biological data also enter into prediction surface
      yrs = 1970:p$current.assessment.year

      for( y in yrs ){
        print(y)
        PS = habitat.db ( DS="complete", year=y, p=p )
        PS = fishing.area.designations( PS, type="planar")
        PS$cfa.factor = NULL
        pH = snowcrab.habitat.db ( DS="habitat", yrs=y, model.type="gam.full", vclass="R0.mass" )
        PS$habitat = pH$fit
        PS$habitat.se = pH$se.fit
        rm (pH ) ; gc()

        ips = NULL
        ips = filter.prediction.locations( DS="default", PS=PS )
        ips = unique( c( 
          ips, 
          filter.prediction.locations( DS="limit.to.near.survey.stations", PS=PS, y=y, p=p ) ) 
        ) # add areas near survey locations even if beyond the range of habitat

        PS = PS[ips, ]

        sa.habitat.total = length(ips )
        temp = tapply( PS$t, INDEX=list(PS$cfa, PS$yr ), FUN=mean, na.rm=T )
        temp.sd = tapply( PS$t, INDEX=list(PS$cfa, PS$yr ), FUN=sd, na.rm=T )
         sa = tapply( PS$t, INDEX=list(PS$cfa, PS$yr ), FUN=length )

        dat = as.data.frame.matrix( cbind( sa, temp, temp.sd) )
        names(dat) = c("sa", "t", "t.sd" )
        dat$cfa = rownames(dat)
       
        all = data.frame( cbind( sa=sa.habitat.total, t=mean(PS$t,na.rm=T), t.sd=sd(PS$t,na.rm=T)) ) 
        all$cfa = "cfaall"

        dat = rbind( dat, all)
        dat$yr = y
        dat$sa = dat$sa * (p$pres*p$pres)

        out = rbind( out, dat )
      }
      PS.hab = out
      save( PS.hab, file=fn, compress=T )

      return  ( PS.hab )
    }
  }



