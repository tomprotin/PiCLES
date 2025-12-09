module DataAssimilation

import PiCLES.custom_structures:DensityObservation
import PiCLES.custom_structures:ExtendedDensityObservation
import Plots as plt
import Statistics as stc

using DataFrames, CSV
using StatsBase, LinearAlgebra

export assimilate!, checkForObservations

function assimilate!(model, Δt)
    nObs = length(model.data_assimilation)

    for iObs in 1:nObs
        AssimTime = model.data_assimilation[iObs].time
        if (model.clock.time+0.001) >= AssimTime && model.data_assimilation[iObs].on
            if typeof(model.data_assimilation[iObs]) == DensityObservation
                nMeshPoints = 1
                meshPoints = [model.data_assimilation[iObs].position_ij[1] model.data_assimilation[iObs].position_ij[2]]
                meshWeights = [1.]
            elseif typeof(model.data_assimilation[iObs]) == ExtendedDensityObservation
                nMeshPoints = length(model.data_assimilation[iObs].positions_ij[:,1])
                meshPoints = model.data_assimilation[iObs].positions_ij
                meshWeights = model.data_assimilation[iObs].mesh_weights
            end
            x = model.data_assimilation[iObs].position_xy[1]
            y = model.data_assimilation[iObs].position_xy[2]
            density = model.data_assimilation[iObs].density
            obsMean = model.data_assimilation[iObs].densityMean
            obsCov = model.data_assimilation[iObs].densityCov


            for meshPoint in 1:nMeshPoints
                i = meshPoints[meshPoint,1]
                j = meshPoints[meshPoint,2]
                locallen = length(model.ParticlesAtNode[i][j])

                if !model.data_assimilation[iObs].evaluateAssimilation
                    if locallen != 0
                        alpha = model.data_assimilation[iObs].evaluateAssimilationAlphaRange*meshWeights[1]
                        evaluate_spectrum(model.ParticlesAtNode[i][j],i,j,[-10,20],[-3,27],101,101,"Epanechnikov",model.plot_savepath,model.clock.time,Δt,"forecast")

                        weights = [(model.ParticlesAtNode[i][j][k][1]) for k in 1:locallen]
                        EnergyTot = sum(weights)
                        weight_f = zeros(locallen)
                        kdePart = zeros(locallen)
                        weight_o_tilde = zeros(locallen)
                        speeds = zeros(locallen,2)
                        centeredSpeeds = zeros(locallen,2)
                        for k in 1:locallen
                            speeds[k,:] = [(model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[2]), (model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[3])]
                        end
                        meanSpeeds = [mean(speeds[:,1],FrequencyWeights(weights)), mean(speeds[:,2],FrequencyWeights(weights))]
                        for k in 1:locallen
                            centeredSpeeds[k,:] = speeds[k,:] - meanSpeeds
                        end
                        centeredSpeedAmplitudes = [sqrt(centeredSpeeds[k,1]^2+centeredSpeeds[k,2]^2) for k in 1:locallen]
                        forecast_variance = var(centeredSpeedAmplitudes,FrequencyWeights(weights),corrected=true)
                        for k in 1:locallen
                            particleVx = model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[2]
                            particleVy = model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[3]
                            weight_f[k] = weights[k] / EnergyTot
                            kdePart[k] = KDE(model.ParticlesAtNode[i][j],particleVx,particleVy,"Epanechnikov",sqrt(forecast_variance))
                            weight_o_tilde[k] = density(particleVx,particleVy)/kdePart[k]
                        end
                        weight_o = weight_o_tilde / sum(weight_o_tilde)
                        pi_tilde = weight_f.^(1/(1+alpha)) .* weight_o.^(alpha/(1+alpha))
                        weight_a =  pi_tilde / sum(pi_tilde) * EnergyTot

                        for k in 1:locallen
                            model.ParticlesAtNode[i][j][k][1] = weight_a[k]
                            # model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[1] = log(weight_a[k])
                        end
                        evaluate_spectrum(model.ParticlesAtNode[i][j],i,j,[-10,20],[-3,27],101,101,"Epanechnikov",model.plot_savepath,model.clock.time,Δt,"corrected")
                        evaluate_spectrum(density,i,j,[-10,20],[-3,27],101,101,model.plot_savepath,model.clock.time,Δt,"observed")
                    end
                else

                    # EVALUATIING THE ASSIMILATION OF THE DATA

                    alphaRange = model.data_assimilation[iObs].evaluateAssimilationAlphaRange
                    nAssim = length(alphaRange)
                    correctedMeans = zeros(nAssim,2)
                    correctedVar = zeros(nAssim)
                    if locallen != 0
                        evaluate_spectrum(model.ParticlesAtNode[i][j],i,j,[-10,20],[-3,27],201,201,"Epanechnikov",model.plot_savepath,model.clock.time,Δt,"0_forecast")
                        evaluate_spectrum(density,i,j,[-10,20],[-3,27],101,101,model.plot_savepath,model.clock.time,Δt,"observed")

                        weights = [(model.ParticlesAtNode[i][j][k][1]) for k in 1:locallen]
                        EnergyTot = sum(weights)
                        speeds = zeros(locallen,2)
                        centeredSpeeds = zeros(locallen,2)
                        for k in 1:locallen
                            speeds[k,:] = [(model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[2]), (model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[3])]
                        end
                        meanSpeeds = [mean(speeds[:,1],FrequencyWeights(weights)), mean(speeds[:,2],FrequencyWeights(weights))]
                        for k in 1:locallen
                            centeredSpeeds[k,:] = speeds[k,:] - meanSpeeds
                        end
                        centeredSpeedAmplitudes = [sqrt(centeredSpeeds[k,1]^2+centeredSpeeds[k,2]^2) for k in 1:locallen]
                        forecast_variance = var(centeredSpeedAmplitudes,FrequencyWeights(weights),corrected=true)

                        for p in 1:nAssim
                            currentAlpha = alphaRange[p]*meshWeights[meshPoint]
                            weight_f = zeros(locallen)
                            kdePart = zeros(locallen)
                            weight_o_tilde = zeros(locallen)
                            for k in 1:locallen
                                particleVx = model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[2]
                                particleVy = model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[3]
                                weight_f[k] = weights[k] / EnergyTot
                                kdePart[k] = KDE(model.ParticlesAtNode[i][j],particleVx,particleVy,"Epanechnikov",sqrt(forecast_variance))
                                weight_o_tilde[k] = density(particleVx,particleVy)/kdePart[k]
                            end
                            weight_o = weight_o_tilde / sum(weight_o_tilde)
                            pi_tilde = weight_f.^(1/(1+currentAlpha)) .* weight_o.^(currentAlpha/(1+currentAlpha))
                            weight_a =  pi_tilde / sum(pi_tilde) * EnergyTot

                            for k in 1:locallen
                                model.ParticlesAtNode[i][j][k][1] = weight_a[k]
                                # model.ParticlesAtNode[i][j][k][3].ODEIntegrator.u[1] = log(weight_a[k])
                            end
                            correctedMeans[p,:] = [mean(speeds[:,1],FrequencyWeights(weight_a)), mean(speeds[:,2],FrequencyWeights(weight_a))]
                            corrected_variance = var(centeredSpeedAmplitudes,FrequencyWeights(weight_a),corrected=true)
                            correctedVar[p] = corrected_variance[1,1]
                            evaluate_spectrum(model.ParticlesAtNode[i][j],i,j,[-10,20],[-3,27],201,201,"Epanechnikov",model.plot_savepath,model.clock.time,Δt,string(p)*"_alpha_"*string(currentAlpha)*"_corrected")
                        end
                        sec=string(Int64(floor((model.clock.time+Δt)/60)))
                        dec=string(Int64(floor(10*((model.clock.time+Δt)/60-floor((model.clock.time+Δt)/60)))))

                        currentPath = model.plot_savepath*"/assimilation/"*sec*","*dec*"/"*string(i)*"_"*string(j)

                        plt.plot(alphaRange,correctedMeans[:,1],size=(1080, 1080),color="blue",
                        title = "analysed speed (on x axis) evolution in alpha",
                        xlabel = "alpha",
                        ylabel = "Vx (m/s)")
                        plt.hline!([meanSpeeds[1],meanSpeeds[1]], linestyle=:dash)
                        plt.hline!([obsMean[1],obsMean[1]], linestyle=:dash)
                        plt.savefig(currentPath*"/Vx_evolution.png")

                        plt.plot(alphaRange,correctedMeans[:,2],size=(1080, 1080),color="blue",
                        title = "analysed speed (on y axis) evolution in alpha",
                        xlabel = "alpha",
                        ylabel = "Vy (m/s)")
                        plt.hline!([meanSpeeds[2],meanSpeeds[2]], linestyle=:dash)
                        plt.hline!([obsMean[2],obsMean[2]], linestyle=:dash)
                        plt.savefig(currentPath*"/Vy_evolution.png")

                        plt.plot(alphaRange,correctedVar[:],size=(1080, 1080),color="blue",
                        title = "analysed speed variance evolution in alpha",
                        xlabel = "alpha",
                        ylabel = "sigma")
                        plt.hline!([forecast_variance,forecast_variance], linestyle=:dash)
                        plt.hline!([obsCov[1,1],obsCov[1,1]], linestyle=:dash)
                        plt.savefig(currentPath*"/sigma_evolution.png")

                    end
                end
            end
            model.data_assimilation[iObs].on = false
        end
    end
end

function checkForObservations(model)
    """
    WIP

    temporary method to create a virtual observation to test the assimilation
    """

    nObs = length(model.data_assimilation)
    if nObs == 0
        t = 18.0*60
        i = 78
        j = 28
        xtemp = i * 22.5e4/151
        ytemp = j * 22.5e4/151

        mean = [4.75, 21.35]
        cov = [3.8 -0.65; -0.65 0.17]
        cov1 = inv(cov)

        alphaRange=0:0.5:2.0
        # alphaRange=[1]

        # newObs = DensityObservation((i,j), (xtemp, ytemp), (vx,vy) -> gaussian2D_PDF(vx,vy,mean,cov,cov1), mean, cov, true, t, true, Vector(alphaRange))
        newObs = ExtendedDensityObservation([i j; i+1 j], [1,0.5], (xtemp, ytemp), (vx,vy) -> gaussian2D_PDF(vx,vy,mean,cov,cov1), mean, cov, true, t, true, Vector(alphaRange))

        append!(model.data_assimilation, [newObs])
    end
end

function gaussian2D_PDF(x,y,mean,cov,cov1)
    return 1/(2*pi*det(cov)^0.5) * exp(-0.5*transpose([x,y]-mean)*cov1*([x,y]-mean))
end

function KDE(particleList,X,Y,kernel,sigma)
    if kernel == "Epanechnikov"
        f = EpanechnikovKernel_approx
    elseif kernel == "Gaussian"
        f= GaussianKernel_approx
    end
    nParticles = length(particleList)
    energies = [particleList[k][1] for k in 1:nParticles]
    particleDistances = zeros(nParticles)
    particleWeights = zeros(nParticles)
    for k in 1:nParticles
        particleDistances[k] = sqrt((particleList[k][3].ODEIntegrator.u[2]-X)^2+(particleList[k][3].ODEIntegrator.u[3]-Y)^2)
        particleWeights[k] = f(particleDistances[k], sigma) * energies[k]
    end

    res = sum(particleWeights) / sum(energies)

    # @info "sum = "
    # @info res

    return res
end

function EpanechnikovKernel_approx(d,sigma)
    if d < sigma
        return (1-(Float64(d)/Float64(sigma))^2)
    else
        return 0
    end
end

function GaussianKernel_approx(d,sigma)
    if d < 6*sigma
        return exp(-0.5*(d/sigma)^2)
    else
        return 0
    end
end

function evaluate_spectrum(particleList::Array,i,j,rangeX,rangeY,resolutionX,resolutionY,kernel,save_path,time,Δt,plotTitle)
    nParticles = length(particleList)
    if nParticles != 0
        xGrid = (0:(resolutionX-1)) ./(resolutionX-1) .*(rangeX[2]-rangeX[1]) .+rangeX[1]
        yGrid = (0:(resolutionY-1)) ./(resolutionY-1) .*(rangeY[2]-rangeY[1]) .+rangeY[1]
        Δx = xGrid[2]-xGrid[1]
        Δy = yGrid[2]-yGrid[1]
        sigma = 1.5*(Δy+Δx)/2
        spectrum = zeros(resolutionX,resolutionY)
        for k in 1:resolutionX
            x = xGrid[k]
            for l in 1:resolutionY
                y = yGrid[l]
                spectrum[k,l] = KDE(particleList,x,y,kernel,sigma)
            end
        end
        spectrum = spectrum ./ sum(spectrum)
        sec=string(Int64(floor((time+Δt)/60)))
        dec=string(Int64(floor(10*((time+Δt)/60-floor((time+Δt)/60)))))

        currentPath = mkpath(save_path*"/assimilation/"*sec*","*dec*"/"*string(i)*"_"*string(j))
        CSV.write(currentPath*"/forecast_spectrum_"*sec*","*dec*"_at_"*string(i)*","*string(j)*".csv",DataFrame(spectrum,:auto))
        CSV.write(currentPath*"/speeds"*sec*","*dec*"_at_"*string(i)*","*string(j)*".csv", DataFrame(x=xGrid,y=yGrid))

        p1 = plt.heatmap(xGrid, yGrid, transpose(spectrum),aspect_ratio=:equal,size=(1080, 1080))
        plt.plot!(legend=:none,
        title = plotTitle*" energy spectrum in wave speed",
        xlabel = "Vx (m/s)",
        ylabel = "Vy (m/s)",
        xlims = (rangeX[1],rangeX[2]),
        ylims = (rangeY[1],rangeY[2]))
        plt.hline!([0 0], color = "red", linestyle=:dash)
        plt.vline!([0 0], color = "red", linestyle=:dash)
        plt.savefig(currentPath*"/"*plotTitle*"_energy_spectrum.png")
    else
        @info "Can't evaluate spectrum : no particles present !"
    end
end

function evaluate_spectrum(PDF::Function,i,j,rangeX,rangeY,resolutionX,resolutionY,save_path,time,Δt,plotTitle)
    xGrid = (0:(resolutionX-1)) ./(resolutionX-1) .*(rangeX[2]-rangeX[1]) .+rangeX[1]
    yGrid = (0:(resolutionY-1)) ./(resolutionY-1) .*(rangeY[2]-rangeY[1]) .+rangeY[1]
    Δx = xGrid[2]-xGrid[1]
    Δy = yGrid[2]-yGrid[1]
    spectrum = zeros(resolutionX,resolutionY)
    for k in 1:resolutionX
        x = xGrid[k]
        for l in 1:resolutionY
            y = yGrid[l]
            spectrum[k,l] = PDF(x,y)
        end
    end
    spectrum = spectrum ./ sum(spectrum)
    sec=string(Int64(floor((time+Δt)/60)))
    dec=string(Int64(floor(10*((time+Δt)/60-floor((time+Δt)/60)))))

    currentPath = mkpath(save_path*"/assimilation/"*sec*","*dec*"/"*string(i)*"_"*string(j))
    CSV.write(currentPath*"/forecast_spectrum_"*sec*","*dec*"_at_"*string(i)*","*string(j)*".csv",DataFrame(spectrum,:auto))
    CSV.write(currentPath*"/speeds"*sec*","*dec*"_at_"*string(i)*","*string(j)*".csv", DataFrame(x=xGrid,y=yGrid))

    p1 = plt.heatmap(xGrid, yGrid, transpose(spectrum),aspect_ratio=:equal,size=(1080, 1080))
    plt.plot!(legend=:none,
    title = plotTitle*" energy spectrum in wave speed",
    xlabel = "Vx (m/s)",
    ylabel = "Vy (m/s)",
    xlims = (rangeX[1],rangeX[2]),
    ylims = (rangeY[1],rangeY[2]))
    plt.hline!([0 0], color = "red", linestyle=:dash)
    plt.vline!([0 0], color = "red", linestyle=:dash)
    plt.savefig(currentPath*"/"*plotTitle*"_energy_spectrum.png")
end

end