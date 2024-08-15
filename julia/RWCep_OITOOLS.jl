using OITOOLS
using PyPlot

oifitsfile="../oifits_data/MYSTIC_L2.2023Aug04.RW_Cep.MIRCX_IDL.2023Dec25.AVG10m.oifits" 
pixsize = 0.05 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)

#Read oifits data
data = readoifits(oifitsfile)[1,1];

#Setup for limb darkening power law fit
ft = setup_nfft(data, nx, pixsize);
weights=[1.0,1.0,1.0]
disc = create_component(type="ldpow", name="disc"); #limb darkening power law fit
disc.vis_params[1].val=2.5
disc.vis_params[1].minval=1.5
disc.vis_params[1].maxval=3.0
model=create_model(disc)

#Actual fit
minf, minx, cvis_model, result = fit_model_ultranest(data, model, weights=weights); #interesting local minimum with v2 only!
minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=weights);
mask=vec(disk(npix=nx, diameter=minx[1]/pixsize+1))
prior=vec(model_to_image(model, nx=nx, pixsize=pixsize)).*mask

regularizers = [["centering", 1e5,1e-4], ["compactness", 4e3,1e-4], ["l1l2", 7e6, 1e-4]];


x = reconstruct(prior.*mask, data, ft, regularizers = regularizers, verb = false, maxiter=500, weights=weights);
x = reconstruct(x.*mask, data, ft, regularizers = regularizers, verb = false, maxiter=500, weights=weights);
chi2 = chi2_nfft_f(vec(x), ft, data, weights=weights, verb=true);

#reshape 1d into 2d so that we can display 
y=reshape(x, nx, nx)
imdisp(x, pixsize=pixsize)

#plot residuals
v2_model, t3amp_model, t3phi_model = image_to_obs(x, ft, data)
plot_v2_residuals(data, image_to_v2(x, data, ft), logplot=true)
plot_t3phi_residuals(data, t3phi_model)

#plot uv 
uvplot(data, color="wavelength");