;;; Matt's IDL code to solve VPH grating efficiency, resolution


common grating,thick,rho,nref,deltan,phi,nglass,nair

pro setGrating,thick0,rh0,nref0,deltan0,tilt0,nglass0
common grating
print,"Remember:"
print,"thickness of grating in um"
print,"rho, number of grooves is in  in lines/um"
print,"nref, the index of refraction of the gelatin"
print,"deltan, is the variation in the index of refraction of the gelatin"
print,"phi in degree tilt from the normal"
print,"glass index of refraction, (usually not used)."
thick = thick0
rho = rh0
nref = nref0
deltan = deltan0
phi = (90.+tilt0)/180.*3.14159
nglass = nglass0
nair = 1.0
end

function K
common grating
return,2*3.14159*rho
end

function beta,lambda
common grating
return,2*3.14159*nref/lambda
end

function cs,lambda,theta
common grating
return,cos(theta)-K()/beta(lambda)*cos(phi)
end

function cr,theta
common grating
return,cos(theta)
end

function nus,lambda,theta
common grating
return, 3.14159 * deltan * thick / lambda / sqrt(cs(lambda,theta)*cr(theta))
end

function nup,lambda,theta
common grating
return,nus(lambda,theta)*(-cos(2*(theta-phi)))
end

function xi,lambda,theta
common grating
return, thick*(K()*cos(theta-phi)-K()*K()*lambda/(4*3.14159*nref))/(2*cs(lambda,theta))
end

function etas,lambda,theta
common grating
xi_ = xi(lambda,theta)
nus_ = nus(lambda,theta)
return, (sin(sqrt(nus_^2 + xi_^2)))^2/(1+xi_^2/nus_^2)
end

function etap,lambda,theta
common grating
xi_ = xi(lambda,theta)
nup_ = nup(lambda,theta)
return, (sin(sqrt(nup_^2 + xi_^2)))^2/(1+xi_^2/nup_^2)
end

function eta,lambda,theta
common grating
return, 0.5*(etap(lambda,theta)+etas(lambda,theta)) 
end

function gratangle,angle
common grating
return, asin(sin(angle/180. * 3.14159)/nref)
end


function bragg,lambda
common grating
return, asin(lambda*rho/2.)
end


function rad2deg,angle
return, (angle/3.14159)*180.
end

function deg2rad,angle
 return, (angle/180.)*3.14159
end


function frs, th
common grating
thi = deg2rad(th)
;if abs(nref*thi) gt 1 then return,0
tht = asin( sin(thi)/nref)
return, 1-(sin(tht-thi)/sin(tht+thi))^2
end

function frp, th
common grating
thi = deg2rad(th)
tht = asin(sin(thi)/nref)
return, 1-(tan(tht-thi)/tan(tht+thi))^2
end

function fr,thi
return,0.5*(frs(thi)+frp(thi))
end

function outangle,lmb,th
common grating
blar = rad2deg(asin(lmb * rho - sin(deg2rad(th))))
return, blar
end


function fresnel_T_s,thi,tho
return, 1-((sin(thi-tho)/sin(thi+tho))^2)
end

function fresnel_T_p,thi,tho
return, 1-((tan(thi-tho)/tan(thi+tho))^2)
end

function fresnel_T,thi,tho
return, 0.5 * ( fresnel_T_s(thi,tho) + fresnel_T_p(thi,tho))
end

function snells_law,thi,ni,no
return, asin(sin(thi)*ni/no)
end

function grating_equation,lmb,thi,print=print
common grating
;return, asin(lmb * rho - sin(thi))
delta = phi-3.14159/2. 
th = gratangle(rad2deg(thi))
if keyword_set(print) then begin
print,"phi",phi
print,"thi",th
endif
return, asin(nref * sin( asin(lmb*rho/nref - sin(th + delta)) + delta))
end 

function bragg_angle,lmb
common grating
delta = phi-3.14159/2. 
;print,cos(delta),lmb, rho, rad2deg(asin(lmb*rho/cos(delta))),lmb*rho/2./cos(delta)
;return, asin (lmb * rho /2./ cos(delta))
return, asin ( nref * sin ( (asin(lmb * rho/2./nref)+delta)))
end


function transmission,lmb,thi, no_coating = no_coating, fixed=fixed
common grating
t = 1.0
thia = deg2rad(thi)
thig = snells_law(thia,nair,nglass)
thir = snells_law(thia,nair,nref)
thoa = grating_equation(lmb,thia)
thog = snells_law(thoa,nair,nglass)
thor = snells_law(thia,nair,nref)


t = t * fresnel_t(thig,thir)
t = t * fresnel_t(thor,thog)

if keyword_set(no_coating)  then begin
t = t * fresnel_t(thog, thoa)
t = t * fresnel_t(thia,thig)
endif

if keyword_set(fixed) then t = t*0.95

return,t

end



; BEGIN MAIN

; start the code!


;; set which grating you want to use here:
  use_grating=1
  gtg=use_grating


pdf = 1.0
physics = 0

skylines = [0.3747, 0.4103, 0.4342, 0.4863, 0.4960, 0.5008, 0.3968, 0.3933]
skylinename = ['OIII', 'H_d', 'H_g' , 'H_b', 'OIII', 'OIII', 'CaH', 'CaK']

ballmer0 = [375.0,395.0]
ballmer1 = [405.0, 425.0]

   _lambdamin = 0.500d
   _lambdamax = 1.100d
   _lambda = (indgen(800)/800.0d)*(_lambdamax-_lambdamin)+_lambdamin
; -------------------- Grating prescriptions
; KCRM RL -- 01 May 2018
; Grating #1
if gtg eq 1 then begin
   _name = "KCRM_RL_v1"
;   _rho = 0.514d
   _rho = 0.5d
   _thick = 7.07d
   _nref = 1.35d
;   _deltan = 0.044d*5.0d/_thick
   _deltan = 0.05
   _tilt = 2.2d
;   _tilt = 1.5d
   _nglass = 1.51d  ; BK7 index of refraction
   _nair = 1.0d
   _angles = [11.0d, 14.0d, 16.0d, 19.5d];[9.37d, 12.3d, 14.4d, 15.0d, 16.0d, 18.0d, 19.5d]
   _waveband = [0.530d, 1.08d]
   _dplotband = [-0.015,+0.015]
   _plotband = _waveband+_dplotband
endif

; KCRM RM1 -- 07 May 2018
; Grating #2
if gtg eq 2 then begin
  _name = "KCRM_RM1_v1"
  _rho = 1.22d
  _thick = 8.42d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.04
  _tilt = 2.38d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [22.0d, 25.0d, 28.5d, 32.0d, 34.5d];, 35.0d];[22.0d, 24.0d, 27.5d, 29.5d, 33.0d, 35.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif

; KCRM RM2 -- 07 May 2018
; Grating #3
if gtg eq 3 then begin
  _name = "KCRM_RM2_v1"
  _rho = 0.921d
  _thick = 9.47d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.05
  _tilt = 2.38d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [22.0d, 25.0d, 28.5d, 32.0d, 34.5d];[22.0d, 24.0d, 27.5d, 29.5d, 33.0d, 35.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif

; KCRM RH1 -- 07 May 2018
; Grating #4
if gtg eq 4 then begin
  _name = "KCRM_RH1_v1"
  _rho = 2.420d
  _thick = 8.88d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.08
  _tilt = 1.0d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [42.0d, 45.0d, 50.0d, 53.0d, 56.0d];[38.0d, 42.0d, 45.0d, 46.0d, 50.0d, 53.5d, 55.0d, 56.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif

; KCRM RH2 -- 07 May 2018
; Grating #5
if gtg eq 5 then begin
  _name = "KCRM_RH2_v1"
  _rho = 2.030d
  _thick = 8.98d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.09
  _tilt = 1.0d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [42.0d, 45.0d, 50.0d, 53.0d, 56.0d];[38.0d, 42.0d, 45.0d, 46.0d, 50.0d, 53.5d, 55.0d, 56.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif

; KCRM RH3 -- 07 May 2018
; Grating #6
if gtg eq 6 then begin
  _name = "KCRM_RH3_v1"
  _rho = 1.705d
  _thick = 10.48d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.09
  _tilt = 1.85d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [42.0d, 45.0d, 50.0d, 53.0d, 56.0d];[38.0d, 42.0d, 45.0d, 46.0d, 50.0d, 53.5d, 55.0d, 56.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif

; KCRM RH4 -- 07 May 2018****** NEEDS VETTING
; Grating #7
if gtg eq 7 then begin
  _name = "KCRM_RH4_v2"
  _rho = 1.495d
;  _rho = 1.435d
  _thick = 13.61d
  _nref = 1.45d
  ;   _deltan = 0.044d*5.0d/_thick
  _deltan = 0.09
  _tilt = 1.8d
  _nglass = 1.51d  ; BK7 index of refraction
  _nair = 1.0d
  _angles = [42.0d, 45.0d, 50.0d, 53.0d, 56.0d]
  _waveband = [0.530d, 1.08d]
  _dplotband = [-0.015,+0.015]
  _plotband = _waveband+_dplotband
endif



print,"Hello"
print,_name
print,gtg
   set_plot,'WIN'

;setgrating,3.6, 3.8, 1.36, 0.14, -1.,1.5


setgrating,_thick,_rho,_nref,_deltan,_tilt,_nglass
print,"yes"

if pdf then begin 
   set_plot,'ps'
   !p.color = fsc_color('black')
   !p.charsize=1.0
   !p.font = 2
   !p.thick = 1.5
   th0 = 5.0
   th1 = 15.0
   device,/color,/encapsulated,/inches,xsize=6.5,ysize=3.75,filename=_name+'.eps'

   textdata=_name+'.fits'


   
   
endif else begin
   !p.background = fsc_color('beige')
   !p.color = fsc_color('black')
   !p.thick=2.0
   !p.charsize=1.5
   th0 = 2.0
   th1 = 7.0
   window,0,xsize=640,ysize=480,xpos=10,ypos=800
endelse

plot,[0,0],[1,1],xrange=1000*_plotband,/xs,xtitle='Wavelength (nm)', ytitle='Theoretical Efficiency', title=_name, yrange=[0,1.15],/ys,pos=[0.1,0.3,0.98,0.93]

code_deg = "260B
;"
letter_deg = "!9"+string(code_deg)+"!X"

code_alpha = "141B
;"
letter_alpha = "!9"+string(code_alpha)+"!X"

code_beta = "142B
;"
letter_beta = "!9"+string(code_beta)+"!X"


x = [ballmer0[0]>_plotband[0]*1000,ballmer0[0]>_plotband[0]*1000,ballmer1[1]<_plotband[1]*1000,ballmer1[1]<_plotband[1]*1000]
y = [!y.crange[0],!y.crange[1],!y.crange[1],!y.crange[0]]
if physics then polyfill,x,y,/data,color=fsc_color("khaki")

createflag=1
trace=fltarr(6,n_elements(_lambda))

for c_ang = 0, n_elements(_angles)-1 do begin
   angle = _angles[c_ang]
   res = eta(_lambda,gratangle(_angles[c_ang]))
   res_p = etap(_lambda,gratangle(_angles[c_ang]))
   res_s = etas(_lambda,gratangle(_angles[c_ang]))
   beta = asin(_lambda*_rho-sin(angle/!radeg))*!radeg
   if c_ang eq 0 then plot,_lambda*1000.0,res,xrange=1000*_plotband,/xs, yrange=[0,1.15],/ys,/noerase,color=fsc_color('black'),pos=[0.1,0.3,0.98,0.93];,title='WHAAAAT'
 ;  oplot,_lambda*1000.0,res,color=fsc_color('black'),thick=5.0
   mx = max(res,mi)
   trace[0,*]=_lambda*1000.0
   trace[1,*]=angle
   trace[2,*]=beta
   trace[3,*]=res
   trace[4,*]=res_p
   trace[5,*]=res_s
   mwrfits,trace,textdata,create=createflag
   createflag=0

   
   for j=30,100 do begin
      vl = mx/100.0*j
      qv = where(res gt vl,nq)
     ; stop
      mnmx = minmax(beta[qv])
      if mnmx[1]-mnmx[0] lt 11.6 then break
   endfor; j
   betamax = beta[mi]
;   q = where(beta gt betamax-5.8 and beta lt betamax+5.8)
   oplot,_lambda[qv]*1000.0,res[qv],color=fsc_color('red'),thick=th1

                                ; let's compute the minimum and
                                ; maximum wavelengths that
                                ; would be put on the CCD here.

   ccd_wvmin = _lambda[min(qv)]
   ccd_wvmax = _lambda[max(qv)]

   ccd_betamin = beta[min(qv)]
   ccd_betamax = beta[max(qv)]
   
   ccd_betacenter = (ccd_betamin + ccd_betamax)/2.0

   ccd_wvcenter = (sin(angle/!radeg)+sin(ccd_betacenter/!radeg))/_rho

   q = where(beta gt betamax-1.93 and beta lt betamax+1.93)
 ;  oplot,_lambda[q]*1000.0,res[q],color=fsc_color('blue'),thick=th1
   oplot,_lambda*1000,res,color=fsc_color('forestgreen'),thick=th0
;   xyouts,/data,ccd_wvcenter*1000.0,res[mi[0]]+0.10,letter_alpha+" = "+string(angle,"(f4.1)")+letter_deg,alignment=0.5,/noclip,charsize=0.75
;   xyouts,/data,ccd_wvcenter*1000.0,res[mi[0]]+0.05,letter_beta+" = "+string(ccd_betacenter,"(f4.1)")+letter_deg,alignment=0.5,/noclip,charsize=0.75
   ;print,14.*150*_lambda[mi[0]]*_rho/cos(angle/!radeg)


   xyouts,/data,_lambda[mi[0]]*1000.0,res[mi[0]]+0.10,letter_alpha+" = "+string(angle,"(f5.2)")+letter_deg,alignment=0.5,/noclip,charsize=0.75
   xyouts,/data,_lambda[mi[0]]*1000.0,res[mi[0]]+0.05,letter_beta+" = "+string(betamax,"(f5.2)")+letter_deg,alignment=0.5,/noclip,charsize=0.75
   print,'Central wave = ',strtrim(_lambda[mi[0]]*1000,1),' nm'
   print,'R(',strtrim(angle,1),'): ',strtrim(14.*150*_lambda[mi[0]]*_rho/cos(angle/!radeg),1)
   



   oplot,[_lambda[mi[0]],_lambda[mi[0]]]*1000.0,[0,res[mi[0]]],linestyle=1,color=fsc_color('purple'),thick=2.0

;   oplot,[ccd_wvcenter,ccd_wvcenter]*1000.0,[0,res[mi[0]]],linestyle=2,color=fsc_color('purple'),thick=2.0




;   print,ccd_wvcenter

endfor 

!p.color = fsc_color('blue')
xrg = !x.crange[1]-!x.crange[0]
x0 = !x.crange[0]

code_rho = "162B
;" 
letter_rho = "!9"+string(code_rho)+"!x"

code_phi = "152B
;" 
letter_phi = "!9"+string(code_phi)+"!x"

code_mu = "155B
;" 
letter_mu = "!9"+string(code_mu)+"!x"

code_lambda = "154B
;" 
letter_lambda = "!9"+string(code_lambda)+"!x"



xyouts,/noclip,/data,x0,-0.26, letter_rho+" = "+string(_rho*1000,"(i4)")+' lines/mm'
xyouts,/noclip,/data,x0+0.0*xrg,-0.46, letter_phi+" = "+string(_tilt,"(f5.1)")+letter_deg
xyouts,/noclip,/data,x0+0.0*xrg,-0.36,"d = "+string(_thick,"(f4.1)")+" "+letter_mu+"m"
xyouts,/noclip,/data,x0+0.4*xrg,-0.26, " n = "+string(nref,"(f4.2)")
xyouts,/noclip,/data,x0+0.4*xrg,-0.36, "dn = "+string(_deltan,"(f4.2)")
xyouts,/noclip,/data,x0+0.75*xrg,-0.26, letter_lambda+"!Dmin!N = "+string(1000*_waveband[0],"(f5.1)")+" nm"
xyouts,/noclip,/data,x0+0.75*xrg,-0.36, letter_lambda+"!Dmax!N = "+string(1000*_waveband[1],"(f6.1)")+" nm"


oplot,1000*[_waveband[0],_waveband[0]],!y.crange,color=fsc_color('red'),linestyle=2,thick=3.0
oplot,1000*[_waveband[1],_waveband[1]],!y.crange,color=fsc_color('red'),linestyle=2,thick=3.0




for k=0,n_elements(skylines)-1 do begin
if physics then    oplot,1000*[skylines[k],skylines[k]],!y.crange,color=fsc_color('blue'),linestyle=1,thick=1.0
endfor; k 

;oplot,[350,350],!y.crange,color=fsc_color('blue'),thick=1.0
;oplot,[560,560],!y.crange,color=fsc_color('blue'),thick=1.0

if pdf then begin
device,/close
set_plot,'WIN'

   ans=mrdfits(textdata,0,hdr)
   sxaddpar,hdr,'NAME',_name,'Grating name'
   sxaddpar,hdr,'RHO',_rho,'Fringe frequency lines/um'
   sxaddpar,hdr,'THICK',_thick,'Gelatin thickness'
   sxaddpar,hdr,'nGEL',_nref,'Gelatin index of refraction'
   sxaddpar,hdr,'dn',_deltan,'Gelatin index of refraction modulation'
   sxaddpar,hdr,'TILT',_tilt,'Fringe tilt'
   sxaddpar,hdr,'nGLASS',_nglass,'Substrate index of refraction'
   sxaddpar,hdr,'nANGLES',n_elements(_angles),'Number of angles in file'
   sxaddpar,hdr,'WAVEMIN',_waveband[0],'Min wavelength in plots'
   sxaddpar,hdr,'WAVEMAX',_waveband[1],'Max wavelength in plots'
   sxaddpar,hdr,'COLUMNS','WAVE, ALPHA, BETA, ETA, ETAP, ETAS'
   modfits,textdata,0,hdr


;spawn,'epstopdf '+_name+'.eps'
endif else begin

endelse


print,"Done"
end 
