function gaussdis(x0,σ,n) #Función para generar un conjunto de n puntos con una distribución gausiana centrada en x0

      # Inicializamos el vector x y definimos un primer valor
      x = zeros(n)
      xc = x0
      xp = 0.0
      #acep = 0.0

      # Iniciamos M-H
      for i = 1:n
            for j = 1:50
                  ϵ = 2.0*(rand(Float64)-0.5)
                  xp = xc+ϵ*2.0*σ
                  pacep=exp(-(xp-x0)^2/(2.0*σ^2.0))/exp(-(xc-x0)^2/(2.0*σ^2.0))
                  r = rand(Float64)
                  if r < pacep
                        xc = xp
                  end
            end
            x[i] = xc
      end
      
      return x
      
end

#=histogram(x,bins = 0:pi/11:2.0*pi)=#



#=function energ(r,l0,Ω,tet0,d) # Programa para calcular la energía de una configuración r con parámetros l_0, Ω ,θ_0

      # Inicializamos la energía y calculamos la matrix de ángulos
      e = 0.0
      #θ = ang(r,d)

      # Calculamos la energía del primer punto
      l = sqrt((r[2,1]-r[1,1])^2+(r[2,2]-r[1,2])^2)
      e += ((l-l0)^2.0)/2.0 
      
      
      for i = 2:d-1
            lp = sqrt((r[i+1,1]-r[i,1])^2.0+(r[i+1,2]-r[i,2])^2.0)
            lm = sqrt((r[i-1,1]-r[i,1])^2.0+(r[i-1,2]-r[i,2])^2.0)
            lb2 = ((r[i+1,1]-r[i,1])-(r[i,1]-r[i-1,1]))^2.0+((r[i+1,2]-r[i,2])-(r[i,2]-r[i-1,2]))^2.0
            e += ((lp-l0)^2.0)/2.0 + ((lm-l0)^2.0)/2.0 + Ω*lb2/(4.0*l0^2.0)  #(Ω/2.0)*(θ[i]-tet0)^2.0
      end

      # Calculamos la energía del último punto
      l = sqrt((r[d,1]-r[d-1,1])^2+(r[d,2]-r[d-1,2])^2)
      e += ((l-l0)^2.0)/2.0 


      return e
end=#



function energ(r,Ω,d,t,τ,elong) # Programa para calcular la energía de una configuración r con parámetros l_0, Ω ,θ_0

      # Inicializamos la energía y calculamos la matrix de ángulos
      e = 0.0
      if elong==true
            l0 = 2.0^((t-1)/τ)
      elseif elong==false
            l0 = 1.0
      else return "ERROR"
      end 
      θ = ang(r,d)
      ϵ = 0.1
      

      # Calculamos la energía del primer punto
      l = sqrt((r[2,1]-r[1,1])^2+(r[2,2]-r[1,2])^2)
      e += ((l-l0)^2.0)/2.0 
      for j = 2:d # Potencial de Lennard-Jones
            l = sqrt((r[j,1]-r[1,1])^2+(r[j,2]-r[1,2])^2)
            e += 4.0*ϵ*((l0/(l*2^(1/6)))^12-(l0/(l*2^(1/6)))^6)
      end
      
      
      for i = 2:d-1
            l = sqrt((r[i+1,1]-r[i,1])^2.0+(r[i+1,2]-r[i,2])^2.0)
            #lb2 = ((r[i+1,1]-r[i,1])-(r[i,1]-r[i-1,1]))^2.0+((r[i+1,2]-r[i,2])-(r[i,2]-r[i-1,2]))^2.0
            e += ((l-l0)^2.0)/2.0 + Ω*(θ[i]-pi)^2.0 #+ Ω*lb2/(4.0*l0^2.0) 
            for j = i:d
                  l = sqrt((r[j,1]-r[i,1])^2+(r[j,2]-r[i,2])^2)
                  e += 4.0*ϵ*((l0/(l*2^(1/6)))^12-(l0/(l*2^(1/6)))^6)
            end
      end

      #e += step(r,d,l0)

      return e 
end



function energver(r,Ω,d,t,τ,elong,j) # Programa para calcular la energía de una configuración r con parámetros l_0, Ω ,θ_0

      # Inicializamos la energía y calculamos la matrix de ángulos
      e = 0.0
      if elong==true
            l0 = 2.0^((t-1)/τ)
      elseif elong==false
            l0 = 1.0
      else return "ERROR"
      end 
      θ = ang(r,d)
      ϵ = 0.1

      if j==1
            l = sqrt((r[2,1]-r[1,1])^2.0+(r[2,2]-r[1,2])^2.0)
            return ((l-l0)^2.0)/2.0
      elseif j==d
            l = sqrt((r[d,1]-r[d-1,1])^2.0+(r[d,2]-r[d-1,2])^2.0)
            return ((l-l0)^2.0)/2.0
      end
            
      
      lp = sqrt((r[j+1,1]-r[j,1])^2.0+(r[j+1,2]-r[j,2])^2.0)
      lm = sqrt((r[j-1,1]-r[j,1])^2.0+(r[j-1,2]-r[j,2])^2.0)
      e = ((lp-l0)^2.0)/2.0 + ((lm-l0)^2.0)/2.0 + Ω*(θ[j]-pi)^2.0  

      #e += step(r,d,l0)

      return e 
end



function enelbe(r,Ω,d,t,τ,elong) # Programa para calcular la energía de una configuración r con parámetros l_0, Ω ,θ_0

      # Inicializamos la energía y calculamos la matrix de ángulos
      enel = 0.0
      enbe = 0.0
      if elong==true
            l0 = 2.0^((t-1)/τ)
      elseif elong==false
            l0 = 1.0
      else return "ERROR"
      end 
      θ = ang(r,d)
      

      # Calculamos la energía del primer punto
      l = sqrt((r[2,1]-r[1,1])^2+(r[2,2]-r[1,2])^2)
      enel += ((l-l0)^2.0)/2.0 
      
      
      for i = 2:d-1
            l = sqrt((r[i+1,1]-r[i,1])^2.0+(r[i+1,2]-r[i,2])^2.0)
            #lb2 = ((r[i+1,1]-r[i,1])-(r[i,1]-r[i-1,1]))^2.0+((r[i+1,2]-r[i,2])-(r[i,2]-r[i-1,2]))^2.0
            enel += ((l-l0)^2.0)/2.0 
            enbe += Ω*(θ[i]-pi)^2.0 #+ Ω*lb2/(4.0*l0^2.0) 
      end

      return enel, enbe
end



function ead(r,e0,ax,ay,r0,d,n)

      for j = 1:d
            eax = 0.0
            eay = 0.0
            for i = 1:n
                  eax += ax[i]*(cos(2.0*pi*i*(r[j,1])/(r0))+sin(2.0*pi*i*(r[j,1])/(r0)))
                  eay += ay[i]*sin(2.0*pi*i*(r[j,2])/(r0))
            end
            e0 += eax + eay
      end

      return   e0   
end

#=
x = range(1, 25, length=100)
y = range(1, 25, length=100)
e = zeros(25,25)
for i = 1:25
      for j = 1:25
            e[i,j]=ead(i,j,0.0,0.2,26,100)
      end
end
=#



function step(r,d,l0)
      e = 0.0
      δ = 0.25
      x0 = d/2.0
      #x1 = d/3.0
      #x2 = d*2.0/3.0
      for i = 1:d
            #=if r[i,1]<x0
                  e += -1
            elseif r[i,1]>=x0
                  e += -2
            end=#
            e += -δ/(1+exp(-(r[i,1]-x0)/l0))
            #e += -1/(1+exp(-(r[i,1]-x1)/l0)) + 1/(1+exp(-(r[i,1]-x2)/l0))
      end

      return e
end



function cevol!(r,β,Ω,d,t,τ,elong)#=,ead0,ax,ay,r0,n=#

      acep = 0.0
      if elong==true
            l0 = 2.0^((t-1)/τ)
      elseif elong==false
            l0 = 1.0
      else return "ERROR"
      end 


      for i = 1:2*d
                  
            # Elegimos una coordenada aleatoria
            j = floor(Int,rand(Float64)*d) + 1

            # Calculamos la energía de la configuración anterior
            e0 = energ(r,Ω,d,t,τ,elong) #+ ead(r,ead0,ax,ay,r0,d,n)

            # Actualizamos la coordenada
            δ = l0*2.0*(rand(Float64)-0.5)/20.0 # Valor aleatorio en [-l0/c,+l0/c]
            σ = 2.0*pi*rand(Float64) # Ángulo aleatorio en [0,2*π)
            r[j,1] += δ*cos(σ)
            r[j,2] += δ*sin(σ)

            # Calculamos la energía de la nueva configuración
            ep = energ(r,Ω,d,t,τ,elong) #+ ead(r,ead0,ax,ay,r0,d,n)

            pacep = exp(-β*(ep-e0))
            pr = rand(Float64)
            #println(pacep,pr)

            # Aplicamos la selección de M-H y comprobamos si la nueva configuración tiene loops cerrados
            loop = loops(r,d)
            if pr > pacep
                  r[j,1] -= (δ)*cos(σ)
                  r[j,2] -= (δ)*sin(σ)
            else acep += 1.0
            end                 
      end

      return acep/d
end



function ang(r,d)

      θ = zeros(d)
      θ[1], θ[d] = pi, pi # No se van a usar pero por precauzión se definen en equilibrio

      for i = 2:d-1
            prod = (r[i+1,1]-r[i,1])*(r[i-1,1]-r[i,1])+(r[i+1,2]-r[i,2])*(r[i-1,2]-r[i,2])  
            norm = sqrt((r[i+1,1]-r[i,1])^2+(r[i+1,2]-r[i,2])^2)*sqrt((r[i-1,1]-r[i,1])^2+(r[i-1,2]-r[i,2])^2)
            
            if  norm == 0.0 
                  θ[i] = 0.0
            else 
                  θ[i] = acos(cor(prod/norm))
            end
      end
      

      return θ      
end



function cor(x)

      if x > 1.0
            return 1.0
      elseif x < -1.0
            return -1.0
      end

      return x
end



function loops(r,d)

      θ = ang(r,d)
      loop = 0.0

      for i = 2:d-1
            loop += phimod(pi-θ[i]) 
            if loop >= 2.0*pi
                  return true
            end
      end

      return false
      
end



function phimod(α)
      α = mod(α,2.0*pi)
      if α > pi
            return α - 2.0*pi
      elseif α < -pi
            return α + 2.0*pi
      end

      return α
end



function curv(r,d)
      c = 0.0
      xp = 0.0
      yp = 0.0
      xpp = 0.0 
      ypp = 0.0
      for i = 2:d-1
            xp = (r[i+1,1] - r[i-1,1])/2.0
            yp = (r[i+1,2] - r[i-1,2])/2.0
            xpp = (r[i+1,1] - 2.0*r[i,1] + r[i-1,1])
            ypp = (r[i+1,2] - 2.0*r[i,2] + r[i-1,2])
            c += abs(xp*ypp + yp*xpp)/(xp^2.0 + yp^2.0)^(1.5)
      end

      c = c/(d-2)

      return c      
end



function bacteri!(Ω,β,d,f)

      a = 0.0
      touch("conf.txt")
      conf = open("conf.txt","w")
      #=ead0 = 0.0
      n = 5
      σ = 0.2
      ax = zeros(n)
      ax = gaussdis(0.0,σ,n)
      ay = zeros(n)
      ay = gaussdis(0.0,σ,n)
      r0 = 1.=#

      # Generamos un primer vector r como una recta con todos los nodos en la posición de equilibrio
      r = zeros(d,2)
      for i = 2:d 
            r[i,1] = (i-1) + (rand(Float64)-0.5)*(1.0/3.0)
            r[i,2] = (rand(Float64)-0.5)*(1.0/3.0)
      end

      for i = 1:f
            # Realizamos un paso de Montecarlo
            a += cevol!(r,β,Ω,d,i,f,true)#=,ead0,ax,ay,r0,n=#

            #=for j = 1:d
                  write(conf,string(round(r[j,1], sigdigits=5)),string(round(r[j,2], sigdigits=5)))
            end=#
            if mod(i*150,f)==0.0
                  writedlm(conf,r)
            end
      end
      close(conf)

      
      println("La tasa de evolución promedio es: ",a/f)

      return 
end
#=
using DelimitedFiles

x = r[:,1]
y = r[:,2]
scatter(x, y)
plot([e,en,eb],label=["Energía total" "Energía elástica" "Energía bending"])
=#




function parambet(Ω,d,f,n)

      # Iniciamos los vectores
      β = zeros(10)
      bed = zeros(f,10)
      l = zeros(f,10)
      e = zeros(f,10)
      r = zeros(d,2)

      for i = 1:10
            β[i] = 10.0^(log10(15.0) + log10(40.0/15.0)*(i-1)/9.0)

            for j = 1:n

                  for i = 2:d 
                        r[i,1] = (i-1) + 2.0*(rand(Float64)-0.5)/3.0
                        r[i,2] = 2.0*(rand(Float64)-0.5)/3.0
                  end

                  for c = 1:d-1
                        l[1,i] += sqrt((r[c+1,1]-r[c,1])^2.0+(r[c+1,2]-r[c,2])^2.0)
                  end

                  bed[1,i] += ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/n

                  e[1,i] += energ(r,Ω,d,1,f,true)/n

                  for k = 2:f
                        cevol!(r,β[i],Ω,d,k,f,true)

                        for c = 1:d-1
                              l[k,i] += sqrt((r[c+1,1]-r[c,1])^2.0+(r[c+1,2]-r[c,2])^2.0)
                        end

                        bed[k,i] += ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/n

                        e[k,i] += energ(r,Ω,d,k,f,true)/n
                  end
            end
      end

      l .= l./((d-1)*n)

      for i = 1:10
            for k = 1:f
                  
                  bed[k,i] = sqrt(bed[k,i])/(d-1)
                                         
            end

      end

      return β, bed, l, e
    
end

#=a1 = "β = "*string(b[1])
a2 = "β = "*string(b[2])
a3 = "β = "*string(b[3])
a4 = "β = "*string(b[4])
a5 = "β = "*string(b[5])
a6 = "β = "*string(b[6])
a7 = "β = "*string(b[7])
a8 = "β = "*string(b[8])
a9 = "β = "*string(b[9])
a10 = "β = "*string(b[10])

m1=plot(msvn,label=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10],title="√<r^2>/L dependence on β")
m2=plot(msvn[:,5:10],label=[a5 a6 a7 a8 a9 a10],title="√<r^2>/L dependence on β",yrange=[0.98,1.02])

e1=plot(e,label=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10],title="Energy decay dependence on β")
e2=plot(e[:,5:10],label=[a5 a6 a7 a8 a9 a10],title="Energy decay dependence on β",yrange=[0,10])
=#



function betau(Ω,d,f,n)

      # Iniciamos los vectores
      p = 10
      β = zeros(p)
      betin = 100.0
      betfin = 500.0
      e = zeros(f,p)
      e_cut = 0.0
      τ = zeros(p)
      r = zeros(d,2)

      for i = 1:p
            β[i] = 10.0^(log10(betin) + log10(betfin/betin)*(i-1)/(p-1))
            conf = 0.0

            for j = 1:n
                  a = 0

                  for i = 1:d 
                        r[i,1] = (i-1) + 2.0*(rand(Float64)-0.5)/10.0
                        r[i,2] = 2.0*(rand(Float64)-0.5)/10.0
                  end
                  e_cut = energ(r,Ω,d,1,f,false)/exp(1)

                  e[1,i] += energ(r,Ω,d,1,f,false)/n

                  for k = 2:f
                        cevol!(r,β[i],Ω,d,k,f,false)
                        ei = energ(r,Ω,d,k,f,false)

                        if a == 0 && ei < e_cut
                              τ[i] += k-0.5
                              a += 1
                        end

                        e[k,i] += ei/n
                  end
                  conf += a
            end
            println(raw"El tiempo característico para β = ",β[i]," es ",τ[i]/n," (",conf/n,", ",τ[i]/conf,")")
      end

      τ .= τ./n

      return β, e, τ
    
end



function paramom(β,d,f,n)
      # Iniciamos los vectores
      p = 5
      Ω = zeros(p)
      msv = zeros(f,p)
      l = zeros(f,p)
      e = zeros(f,p)

      # Generamos un primer vector r como una recta con todos los nodos en la posición de equilibrio y lo perturbamos ligeramente
      r0 = zeros(d,2)
      for i = 2:d 
            r0[i,1] = (i-1) + (rand(Float64)-0.5)*(1.0/3.0)
            r0[i,2] = (rand(Float64)-0.5)*(1.0/3.0)
      end

      r = zeros(d,2)

      println("Iniciando cálculo de los valores de cada configuración")     
      for i = 1:p
            Ω[i] = 10.0^(-2.0 + 2.0*(i-1)/(p-1))
            println("Iniciando cálculo para Ω = ",Ω[i])
            for j = 1:n
                  r .= r0
                  for c = 1:d-1
                        l[1,i] += sqrt((r[c+1,1]-r[c,1])^2.0+(r[c+1,2]-r[c,2])^2.0)
                  end

                  msv[1,i] += ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/n

                  e[1,i] += energ(r,Ω[i],d,1,f,true)/n
                  for k = 2:f
                        cevol!(r,β,Ω[i],d,k,f,true)

                        for c = 1:d-1
                              l[k,i] += sqrt((r[c+1,1]-r[c,1])^2.0+(r[c+1,2]-r[c,2])^2.0)
                        end

                        msv[k,i] += ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/n

                        e[k,i] += energ(r,Ω[i],d,k,f,true)/n
                  end
                  if mod(10*j,n)==0
                        print(" · ")
                  end
            end
            println(" ")
      end
      
      l .= l./((d-1)*n)

      println("Iniciando cálculo de los valores promedio y sus errores")     
      for i = 1:p
            for k = 1:f
                  
                  msv[k,i] = sqrt(msv[k,i])/(d-1)
                                         
            end

      end

      return Ω, msv, l, e
    
end

#=a1 = "Ω = "*string(o[1]);
a2 = "Ω = "*string(o[2]);
a3 = "Ω = "*string(o[3]);
a4 = "Ω = "*string(o[4]);
a5 = "Ω = "*string(o[5]);
a6 = "Ω = "*string(o[6]);
a7 = "Ω = "*string(o[7]);
a8 = "Ω = "*string(o[8]);
a9 = "Ω = "*string(o[9]);
a10 = "Ω = "*string(o[10]);

m1=plot(msvn,label=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10],title="√<r^2>/L dependence on Ω")
m2=plot(msvn[:,5:10],label=[a5 a6 a7 a8 a9 a10],title="√<r^2>/L dependence on Ω",yrange=[0.98,1.02])

e1=plot(e,label=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10],title="Energy decay dependence on Ω")
e2=plot(e[:,5:10],label=[a5 a6 a7 a8 a9 a10],title="Energy decay dependence on Ω",yrange=[0,10])
=#



function li(Ω,β,d,f,n)

      # Iniciamos los vectores
      li = zeros(f,n)
      l = zeros(f)
      sl = zeros(f)

      ri2 = zeros(f,n)
      r2 = zeros(f)
      sr = zeros(f)

      ci = zeros(f,n)
      c = zeros(f)
      sc = zeros(f)

      r = zeros(d,2)

      println("Iniciando cálculo de los valores de cada configuración")     
      for j = 1:n

            for i = 1:d 
                  r[i,1] = (i-1) + (rand(Float64)-0.5)*(1.0/3.0) # x_i
                  r[i,2] = (rand(Float64)-0.5)*(1.0/3.0) # y_i
            end

            for i = 1:d-1
                  li[1,j] += sqrt((r[i+1,1]-r[i,1])^2.0+(r[i+1,2]-r[i,2])^2.0)
            end

            li[1,j] = li[1,j]/(d-1)
            ri2[1,j] = ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/(d-1)^2.0
            ci[1,j] += curv(r,d)
            
            

            for t = 2:f
                  cevol!(r,β,Ω,d,t,f,true)
                  for i = 1:d-1
                        li[t,j] += sqrt((r[i+1,1]-r[i,1])^2.0+(r[i+1,2]-r[i,2])^2.0)
                  end
                  li[t,j] = li[t,j]/(d-1)
                  ri2[t,j] = ((r[d,1]-r[1,1])^2.0+(r[d,2]-r[1,2])^2.0)/(d-1)^2.0
                  ci[t,j] += curv(r,d)
            end
            if mod(10*j,n)==0
                  print(" · ")
            end
      end
      println(" ")

      println("Iniciando cálculo de los valores promedio y sus errores")     

      for t = 1:f
            for i = 1:n
                  l[t] += li[t,i]
                  r2[t] += ri2[t,i]
                  c[t] += ci[t,i]
            end

            l[t] = l[t]/n
            r2[t] = sqrt(r2[t]/n)
            c[t] = c[t]/n

            for i = 1:n
                  sl[t] += (li[t,i] - l[t])^2.0
                  sr[t] += (sqrt(ri2[t,i]) - r2[t])^2.0
                  sc[t] += (ci[t,i] - c[t])^2.0
            end

            sl[t] = sqrt(sl[t])/n
            sr[t] = sqrt(sr[t])/n
            sc[t] = sqrt(sc[t])/n
      end

      return l, sl, r2, sr, c, sc
end
#=plot(c,ribbon=sc,fillalpha=.5)
l = zeros(1650,3);
c = zeros(1650,3);
for i = 1:3
      l[:,i],_,_,_,c[:,i],_ = li(0.025*4^(i-1),400.0/4^(i-1),16*2^(i-1),1650,1000)
end
pl = plot(l,label=["n=13" "n=26" "n=52"],title="Longitud promedio para varias lattices",legend=:topleft)
pc = plot(c,label=["n=13" "n=26" "n=52"],title="Curvatura promedio para varias lattices",legend=:topleft)
sl = zeros(1650,3);
sc = zeros(1650,3);=#





function ei(Ω,β,d,f,n)

      # Decidimos si el sistema se elonga o no
      elong = true

      # Iniciamos los vectores
      ei = zeros(f,n)
      e = zeros(f)
      se = zeros(f)

      eei = zeros(f,n)
      ee = zeros(f)
      see = zeros(f)

      ebi = zeros(f,n)
      eb = zeros(f)
      seb = zeros(f)

      r = zeros(d,2)

      println("Iniciando cálculo de los valores de cada configuración")     
      for j = 1:n

            for i = 1:d 
                  r[i,1] = (i-1) + (rand(Float64)-0.5)*(1.0/3.0) # x_i
                  r[i,2] = (rand(Float64)-0.5)*(1.0/3.0) # y_i
            end


            ei[1,j] = energ(r,Ω,d,1,f,elong)
            eei[1,j],ebi[1,j] = enelbe(r,Ω,d,1,f,elong)


            for t = 2:f
                  cevol!(r,β,Ω,d,t,f,elong)
                  ei[t,j] = energ(r,Ω,d,t,f,elong)
                  eei[t,j],ebi[t,j] = enelbe(r,Ω,d,t,f,elong)
            end
            if mod(10*j,n)==0
                  print(" · ")
            end
      end
      println(" ")

      println("Iniciando cálculo de los valores promedio y sus errores")     

      for t = 1:f
            for i = 1:n
                  e[t] += ei[t,i]
                  ee[t] += eei[t,i]
                  eb[t] += ebi[t,i]
            end

            e[t] = e[t]/n
            ee[t] = ee[t]/n
            eb[t] = eb[t]/n

            for i = 1:n
                  se[t] += (ei[t,i] - e[t])^2.0
                  see[t] += (eei[t,i] - ee[t])^2.0
                  seb[t] += (ebi[t,i] - eb[t])^2.0
            end

            se[t] = sqrt(se[t])/n
            see[t] = sqrt(see[t])/n
            seb[t] = sqrt(seb[t])/n
      end

      return e, se, ee, see, eb, seb
end

#=
e, se, ee, see, eb, seb = ei(Ω,β,d,f,n)
plot([e,ee,eb],label=["Energía total" "Energía elástica" "Energía bending"])
een = zeros(1650)
ebn = zeros(1650)
een .= ee./e
ebn .= eb./e
plot([een,ebn],label=["Energía elástica normalizada" "Energía bending normalizada"])

x = range(1, τ[i], length=τ[i])
y=@. 2^(x/τ[i])
p1 = plot(x,[y,l],label=["Expected" "Obtained"], title="<l> for τ="*string(τ[i]))
savefig(p,"IQ.pdf")
=#