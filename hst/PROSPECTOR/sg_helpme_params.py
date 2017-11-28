# In [1]
import time, sys, os
tbegin = time.time()
import h5py
import numpy as np
import matplotlib.pyplot as plt

# In [2]
import fsps
import sedpy
import prospect
import emcee
import sg_params as params
from corner import quantile

# Here is some stuff for plotting I may already have but will import again
from matplotlib import cm,axes as ax
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
from sg_likelihood import lnlike_spec, lnlike_phot, write_log
# In [3]
from prospect.models import model_setup
import plot_utils as pread
from prospect.io.read_results import results_from
import xlwt
import xlrd
## stuff for later

# IMMA IMPORT STUFF FOR MY EMAIL BITCHES
import smtplib
from string import Template
from email import encoders

MY_ADDRESS = 'scarlinewgottlieb@gmail.com'
MYPASSWORD = '****'

filename = 'errmessage.txt'
startindex = 0
endindex = 40
def read_template(filename):
    with open(filename, 'r', encoding='utf-8') as template_file:
        template_file_content = template_file.read()
    return Template(template_file_content)

message_template = read_template(filename)
def login():
    # set up the SMTP server
    s = smtplib.SMTP(host='smtp.gmail.com', port = 587)
    #s = smtplib.SMTP(host='smtpout.secureserver.net', port = 465)
    s.starttls()
    s.login(MY_ADDRESS,MYPASSWORD)
    return s

# import necessary packages
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase


def error_email(errobjid, errparam):
    s = login()
    filename = 'errmessage.txt'
    message_template = read_template(filename)

    msg = MIMEMultipart()       # create a message
    part = MIMEBase('application', "octet-stream")
    '''f = 'Result.pdf'
    part.set_payload( open(f,"rb").read() )
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="{0}"'.format(os.path.basename(f)))
    msg.attach(part)'''
     # add in the actual person name to the message template
    message = message_template.substitute(OBJID=errobjid, PARAM=errparam)
    print('Emailing error for ' + errobjid + ' and ' + errparam)
    # setup the parameters of the message
    msg['From']=MY_ADDRESS
    msg['To']='sgottlie@bates.edu'
    msg['Subject']="Eeek! Your code broke!"

    # add in the message body
    msg.attach(MIMEText(message, 'plain'))

    # send the message via the server set up earlier.
    s.send_message(msg)
    
    del msg

def succ_message(objid, parm):
    s = login()
    filename = 'sucmessage.txt'
    message_template = read_template(filename)

    msg = MIMEMultipart()       # create a message
    part = MIMEBase('application', "octet-stream")
     # add in the actual person name to the message template
    message = message_template.substitute(OBJID=objid, PARAM=parm)
    print('Emailing success for ' + objid + ' and ' + parm)
    # setup the parameters of the message
    msg['From']=MY_ADDRESS
    msg['To']='sgottlie@bates.edu'
    msg['Subject']="Never fear! It's alive!"

    # add in the message body
    msg.attach(MIMEText(message, 'plain'))

    # send the message via the server set up earlier.
    s.send_message(msg)
    
    del msg


# In [10]
def lnprobfn(theta):
    """Given a parameter vector, a dictionary of observational data 
    and a model object, return the ln of the posterior. 
    This requires that an sps object (and if using spectra 
    and gaussian processes, a GP object) be instantiated.
    """

    #print('lnprobfn loves pizza')
    # Calculate prior probability and exit if not within prior
    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        #print('oh shit prior')
        return -np.infty
        
    # Generate mean model
    t1 = time.time()
    try:
        spec, phot, x = model.mean_model(theta, obs, sps=sps)
    except(ValueError):
        return -np.infty
    d1 = time.time() - t1
    vectors = {}  # This would be used for noise model weight functions

    # Calculate likelihoods
    t2 = time.time()
    lnp_spec = lnlike_spec(spec, obs=obs, spec_noise=spec_noise)
    lnp_phot = lnlike_phot(phot, obs=obs, phot_noise=phot_noise)
    d2 = time.time() - t2
    if verbose:
        write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

    return lnp_prior + lnp_phot + lnp_spec

# In [11]
from prospect.io import write_results

from prospect import fitting
from scipy.optimize import least_squares

# In [16]
from prospect.likelihood import chi_spec, chi_phot

verbose = False # don't output function calls
tnow = str(int(time.time()))+'/'

###
#paramfile = 'demo_mock_params1.py'
paramfile = 'sg_params.py'
#photfile ='demo_photometry1.dat'

wb = xlwt.Workbook()
ws = wb.add_sheet('Param Results')

# In [4]
clargs = {'param_file':paramfile}
run_params = model_setup.get_run_params(argv=paramfile, **clargs)
run_params['outfile'] = tnow
if not os.path.exists(tnow):
    os.makedirs(tnow)
#print(run_params)

# load photometry from sg_flux.dat
obs = params.load_obs(**run_params)
obs['spectrum'] = None
ws.write(0,0, 'Obj ID')

varnames = ['mass','logzsol','tau','tage','dust2']
coltit = [ty + par for par in varnames for ty in ['low ', 'best ', 'high ']]

for i, lab in enumerate(coltit):
    ws.write(0,i+1, lab)


## HERE BE THE BEGINNING OF THE LOOP
# def prospect(galaxy_name):
with PdfPages(tnow+'total.pdf') as pdf:
    for i in range(startindex,endindex):
    #for i in range(80,240):
        ws.write(i+1,0, obs['objid'][i])
        print(obs['objid'][i])
        print()
        #try: 
        
        obsap = [obs['all_maggies'][j][i] for j in range(0,len(obs['all_maggies']))]
        errap = [obs['all_maggies_unc'][j][i] for j in range(0,len(obs['all_maggies_unc']))]
        obs['maggies'] = tuple(obsap)
        obs['maggies_unc'] = errap

        obs['cmask'] = [obs['all_phot_mask'][j][i] for j in range(0,len(obs['all_phot_mask']))]
        obs['phot_mask'] = obs['cmask']
        obs['mcur'] = obsap
        obs['ecur'] = errap
        if run_params['zred'] != obs['z'][i]:
            run_params['zred'] = obs['z'][i]
            sps = model_setup.load_sps(**run_params)
            spec_noise, phot_noise = model_setup.load_gp(**run_params)
            model = model_setup.load_model(**run_params)

            # Stuff for the file this stuff will go to
        outroot = "{0}_{1}".format(run_params['outfile'], obs['objid'][i])
        try:
            hfilename = outroot + '_mcmc.h5'
            hfile = h5py.File(hfilename, "a")
            #print("Writing to file {}".format(hfilename))
            write_results.write_h5_header(hfile, run_params, model)
            write_results.write_obs_to_h5(hfile, obs)
        except:
            hfile = None
        wspec = sps.csp.wavelengths # *restframe* spectral wavelengths
        a = 1.0 + run_params['zred'] # cosmological redshifting
        # photometric effective wavelengths
        wphot = np.array([f.wave_effective for f in obs['filters']])
        # initial parameters based on `init` values of model_params
        initial_theta = model.rectify_theta(model.initial_theta)
        # generate model
        out_init = model.mean_model(initial_theta, obs, sps=sps) 
        mspec_init, mphot_init, mextra_init = out_init


        # establish bounds
        xmin, xmax = np.min(wphot)*0.8, np.max(wphot)/0.8
        temp = np.interp(np.linspace(xmin,xmax,10000), wspec * a, mspec_init)
        ymin, ymax = temp.min()*0.8, temp.max()/0.8
        ymax = 1e-7


        fig = plt.figure(figsize=(25,12))
        plt.loglog(wspec * a, mspec_init, label='Model spectrum', lw=0.7, color='navy', alpha=0.7)
        plt.errorbar(wphot, mphot_init, label='Model photometry', marker='s',markersize=10, alpha=0.8, ls='', lw=3,markerfacecolor='none', markeredgecolor='blue', markeredgewidth=3)
        plt.errorbar(wphot, obsap, yerr=errap, label='Observed photometry', marker='o', markersize=10, alpha=0.8, ls='', lw=3, ecolor='red', markerfacecolor='none', markeredgecolor='red', markeredgewidth=3)

        # plot Filters
        for f in obs['filters']:
            #print(str(f))
            w, t = f.wavelength.copy(), f.transmission.copy()
            while t.max() > 1:
                t /= 10.
            t = 0.1*(ymax-ymin)*t + ymin
            plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

        # prettify
        plt.xlabel('Wavelength [A]')
        plt.ylabel('Flux Density [maggies]')
        plt.xlim([xmin, xmax])
        #plt.ylim([ymin, ymax])
        plt.legend(loc='best', fontsize=20)
        plt.tight_layout()
        plt.title(obs['objid'][i] + ' Plot 1')
        #pdf.savefig()
        plt.close()




        # start minimization
        min_method = 'levenburg-marquardt'
        nmin = 5 # We'll start from 5 places, 4 of which are drawn from the prior
        ts = time.time()
        pinitial = fitting.minimizer_ball(model.initial_theta.copy(), nmin, model)
        guesses = []
        #print('your mom', pinitial)
        ###
        from sg_likelihood import chi_spec, chi_phot
        def chivecfn(theta):
            """A version of lnprobfn that returns the simple uncertainty 
            normalized residual instead of the log-posterior, for use with 
            least-squares optimization methods like Levenburg-Marquardt.
            """
            lnp_prior = model.prior_product(theta)
            if not np.isfinite(lnp_prior):
                return -np.infty

            # Generate mean model
            t1 = time.time()
            try:
                spec, phot, x = model.mean_model(theta, obs, sps=sps)
            except(ValueError):
                return -np.infty
            d1 = time.time() - t1

            chispec = chi_spec(spec, obs)
            chiphot = chi_phot(phot, obs)
            return np.concatenate([chispec, chiphot])

        ###
        for j, pinit in enumerate(pinitial): #loop over initial guesses
            print('checking initial guesses for ' +varnames[j])
            #print(str(j), pinit)
            try:
                res = least_squares(chivecfn, pinit, method='dogbox', x_scale='jac',
                                xtol=1e-16, ftol=1e-16)
            except:
                #error_email(obs['objid'][i], varnames[j])
                print('error with', obs['objid'][i], varnames[j])
                res = None
            guesses.append(res)


        # Calculate chi-square of the results, and choose the best one
        chisq = [np.sum(r.fun**2) for r in guesses]
        best = np.argmin(chisq)
        initial_center = fitting.reinitialize(guesses[best].x, model,
                                edge_trunc=run_params.get('edge_trunc', 0.1))
        initial_prob = None
        pdur = time.time() - ts

        # output results
        #print('done {0} in {1}s'.format(min_method, pdur))
        #print('best {0} guess: {1}'.format(min_method, initial_center))
        #print('best {0} chi-sq: {1}'.format(min_method, chisq[best]))

        # In [17]
        # initial parameters
        theta = model.rectify_theta(initial_center)
        # generate model
        mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

        # establish bounds
        xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
        temp = np.interp(np.linspace(xmin,xmax,10000), wspec*a, mspec)
        ymin, ymax = temp.min()*0.8, temp.max()/0.8
        ymax = 1e-7

        fig = plt.figure(figsize=(25,12))
        # plot Data, models, and old models
        plt.loglog(wspec * a, mspec_init, label='Old model spectrum',
               lw=0.7, color='gray', alpha=0.5)
        plt.errorbar(wphot, mphot_init, label='Old model Photometry', 
                 marker='s', markersize=10, alpha=0.6, ls='', lw=3, 
                 markerfacecolor='none', markeredgecolor='gray', 
                 markeredgewidth=3)
        plt.loglog(wspec * a, mspec, label='Model spectrum', 
               lw=0.7, color='navy', alpha=0.7)
        plt.errorbar(wphot, mphot, label='Model photometry', 
                 marker='s', markersize=10, alpha=0.8, ls='', lw=3,
                 markerfacecolor='none', markeredgecolor='blue', 
                 markeredgewidth=3)
        plt.errorbar(wphot, obs['mcur'], yerr=obs['ecur'],
                 label='Observed photometry', 
                 marker='o', markersize=10, alpha=0.8, ls='', lw=3, 
                 ecolor='red', markerfacecolor='none', markeredgecolor='red', 
                 markeredgewidth=3)

        # plot filter transmission curves
        for f in obs['filters']:
            w, t = f.wavelength.copy(), f.transmission.copy()
            while t.max() > 1:
                t /= 10.
            t = 0.1*(ymax-ymin)*t + ymin
            plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

        # Prettify
        plt.xlabel('Wavelength [A]')
        plt.ylabel('Flux Density [maggies]')
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        plt.legend(loc='best', fontsize=20)
        plt.tight_layout()
        plt.title(obs['objid'][i] + ' Plot 2')
        #pdf.savefig()

        # In [18] --> this is the step where the code breaks
        postkwargs = {} #any keyord arguments to the lnpostfn would go here.


        tstart = time.time()  # time it

        out = fitting.run_emcee_sampler(lnprobfn, initial_center, model, postkwargs=postkwargs, initial_prob=initial_prob, pool=None, hdf5=hfile, **run_params)

        esampler, burn_p0, burn_prob0 = out
        edur = time.time() - tstart

        #sys.stdout = fout

        #print('done emcee in {0}s'.format(edur))

        # In [19]
        write_results.write_pickles(run_params, model, obs, esampler, guesses,
                                    outroot=outroot, toptimize=pdur, tsample=edur,
                                    sampling_initial_center=initial_center,
                                    post_burnin_center=burn_p0,
                                    post_burnin_prob=burn_prob0)

        if hfile is None:
            hfile = hfilename
        write_results.write_hdf5(hfile, run_params, model, obs, esampler, 
                                 guesses,
                                 toptimize=pdur, tsample=edur,
                                 sampling_initial_center=initial_center,
                                 post_burnin_center=burn_p0,
                                 post_burnin_prob=burn_prob0)

        #print('Finished')

        # grab results, powell results, and our corresponding models
        res, pr, mod = results_from("{}_mcmc.h5".format(outroot))


        # In [22]
        # To see how our MCMC samples look, we can examine a few traces
        choice = np.random.choice
        tracefig = pread.param_evol(res, figsize=(20,10), 
                                    chains=choice(128, size=10, replace=False))
        plt.suptitle(obs['objid'][i] + ' Evolved Params')
        pdf.savefig()

        # In [23]
        # Show samples in a triangle or corner plot
        theta_truth = np.array([run_params[i] 
                                for i in ['mass','logzsol','tau','tage','dust2']])
        theta_truth[0] = np.log10(theta_truth[0])


        fc = res['chain'][:,0::5,:]
        midx = [l=='mass' for l in varnames]
        #fc[:,midx] = np.log10(fc[:,midx])
        fc = np.matrix.transpose(fc)
        print(obs['objid'][i])

        for p00p, things in enumerate(fc):
            q_16, q_50, q_84 = quantile(things, [0.16, 0.5, 0.84],)
            q_m, q_p = q_50-q_16, q_84-q_50

            best = q_50
            low = q_m
            high = q_p
            ws.write(i+1, p00p*3+1, low)
            ws.write(i+1, p00p*3+2, best)
            ws.write(i+1, p00p*3+3, high)


            #print(varnames[p00p], low, best, high)


        cornerfig = pread.subtriangle(res, outname = tnow+obs['objid'][i], start=0, thin=5, truths=theta_truth, fig=plt.subplots(5,5,figsize=(27,27))[0])
        plt.suptitle(obs['objid'][i] + ' Triangle Plot')
        pdf.savefig()

        # In [24]

        # randomly chosen parameters from chain
        randint = np.random.randint
        nwalkers, niter = run_params['nwalkers'], run_params['niter']
        theta = res['chain'][randint(nwalkers), randint(niter)]
        # generate model
        mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)

        # establish bounds
        xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
        temp = np.interp(np.linspace(xmin,xmax,10000), wspec * a, mspec)
        ymin, ymax = temp.min()*0.8, temp.max()/0.8
        fig = plt.figure(figsize=(25,12))

        #fig = plt.figure()

        # plot data and model
        plt.loglog(wspec * a, mspec, label='Model spectrum',
               lw=0.7, color='navy', alpha=0.7)
        plt.errorbar(wphot, mphot, label='Model photometry',
                 marker='s', markersize=10, alpha=0.8, ls='', lw=3, 
                 markerfacecolor='none', markeredgecolor='blue', 
                 markeredgewidth=3)
        plt.errorbar(wphot, obsap, yerr=errap, 
                 label='Observed photometry', ecolor='red', 
             marker='o', markersize=10, ls='', lw=3, alpha=0.8, 
             markerfacecolor='none', markeredgecolor='red', 
             markeredgewidth=3)

    # plot transmission curves
    for f in obs['filters']:
        w, t = f.wavelength.copy(), f.transmission.copy()
        while t.max() > 1:
            t /= 10.
        t = 0.1*(ymax-ymin)*t + ymin
        plt.loglog(w, t, lw=3, color='gray', alpha=0.7)

    plt.xlabel('Wavelength [A]')
    plt.ylabel('Flux Density [maggies]')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.legend(loc='best', fontsize=20)
    plt.title(obs['objid'][i] + ' Final Plot')
    plt.tight_layout()

    pdf.savefig()
    #plt.savefig('test1.png')
    plt.close()
    #plt.savefig('ha.png')
    succ_message(obs['objid'], 'everything')
    '''
        except:
            try:
                error_email(obs['objid'][i], varnames[j])
            except:
                print('You don\'t want an email? FINE. YOU DON\'T GET AN EMAIL, JERK')
                print(obs['objid'][i], varnames[j])
    '''
        
        
wb.save(tnow+'results.xls')

print('HOLY FUCK THAT TOOK FOREVER')

elapsed = (time.time()-tbegin)/3600

print('it took ' + str(elapsed) +' hours!')


