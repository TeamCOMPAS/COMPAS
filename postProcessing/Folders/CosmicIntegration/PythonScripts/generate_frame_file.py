import bilby
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from gwpy.timeseries import TimeSeries, TimeSeriesDict
import random

def multiple_injections (path="/Users/ilyam/Work/COMPASresults/popsynth/Arash/",
                         filename="mergers.txt", dz=0.001, Tobs=1./365.25/24/60, T0=1234567):
    random.seed()
    #path="./"
    input=open(path+filename, 'r')
    input.readline()
    input.readline()
    count=0
    for line in input:
        m1,m2,z,rate=np.float_(line.strip().split("\t"))
        distance = cosmo.luminosity_distance(z).value
        dVcdz = cosmo.differential_comoving_volume(z).value*4*np.pi
        prob_injection = rate * (dVcdz * dz)/1e9 * (Tobs/(1+z))
        assert(prob_injection<1)
        if(random.random()<prob_injection):
            one_injection(m1,m2,z,distance, T0, Tobs)
            count+=count
    input.close()
    print(count)
    return count


def one_injection (m1, m2, z, distance, T0, Tobs):
    print(f'Injecting m1:{m1}, m2:{m2}, z:{z}')
    time=T0+Tobs*365.25*86400*random.random()
    ra=2*np.pi*random.random()
    dec=np.arcsin(-1+2*random.random())
    psi=random.random()*np.pi
    phase=random.random()*2*np.pi
    theta_jn=np.arccos(-1+2*random.random())
    injection_parameters = dict (mass_1=m1*(1+z), mass_2=m2*(1+z), luminosity_distance=distance,
                                 theta_jn=theta_jn, psi=psi, phase=phase, geocent_time=time,
                                 ra=ra, dec=dec, chi_1=0, chi_2=0)
    duration = 32
    sampling_frequency = 2048
    start_time=time-duration
    print(start_time)
    waveform_arguments = dict(waveform_approximant='IMRPhenomPv2', reference_frequency=50, minimum_frequency=40)

# Create the waveform_generator. This is something used for injecting a signal and inference.
    waveform_generator = bilby.gw.WaveformGenerator(
        duration=duration, sampling_frequency=sampling_frequency,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
                                                    waveform_arguments=waveform_arguments)

# Set up two other detectors at Hanford and Livingston. These will use the design sensitivity PSD by default
    interferometers = bilby.gw.detector.InterferometerList(['H1', 'L1'])

# Inject a signal into the network of detectors
    interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency, duration=duration, start_time=start_time)
    interferometers.inject_signal(parameters=injection_parameters,
                              waveform_generator=waveform_generator)


    for interferometer in interferometers:
        signal = interferometer.get_detector_response(waveform_generator.frequency_domain_strain(), injection_parameters)
        #interferometer.plot_data(signal=signal, outdir=path, label='DCO)


#Use gwpy to generate frame file for the network
    H1 = TimeSeries(interferometers[0].time_domain_strain,
                sample_rate=sampling_frequency, unit='strain',
                channel='H1', name='H1')
    L1 = TimeSeries(interferometers[1].time_domain_strain,
                sample_rate=sampling_frequency, unit='strain',
                channel='L1', name='L1')

    ifos = TimeSeriesDict()
    ifos.update(H1=H1, L1=L1)
    ifos.write('frame.gwf')

    return start_time
