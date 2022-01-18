'''
Created on 9/02/2021

@author: macpro
'''
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.client import Client, LSFCluster,SLURMCluster, Future
#from seismicUtilities.AcousticShot import acoustic_shots
from laTest import acousticShot
from utils import parse_args
import time

def main():
    sysargs = parse_args()

    acq = EQUISPACEDAcquisition(boundary=[[0,2000],[0,2000],[0,2000]],
                                dt=0.002,
                                velocity_model=1500,
                                start_source_pos    = [501, 1001],
                                end_source_pos      = [1501, 1001],
                                start_receivers_pos = [[21, 1001]],
                                end_receivers_pos   = [[1981, 1001]],
                                number_of_sources   = 4,
                                number_of_receivers = 4,
                                source_depth = 1899,
                                receivers_depth = 1949)


    acq.add_xml(sysargs.xml)
    acq.limitedAperture(500)
    #acq.calculDt()

    acqs = acq.split(4)

    cluster = SLURMCluster(job_name="seismicAcquisition", nodes=1, cores=4,
                           python="/g/g90/hamon1/geosx/GEOSX/build-quartz-gcc@8.1.0-release/lib/PYGEOSX/bin/python",
                           extra=["-p pdebug"])

    client = Client(cluster)
    client.scale(3)


    maxTime = 0.5
    outputSeismoTraceInterval = 5
    outputWaveFieldInterval = 100

    
    args = []
    for acq in acqs:
        args.append((maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acq,))
        
    #args = (maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acq,)

    
    futures = client.map(acousticShot, args,
                         cmd_line = ["-i", sysargs.xml],
                         cores=2, x_partition=2)
    
    #future = client.submit(acousticShot, args,
    #                       cmd_line=["-i", sysargs.xml],
    #                       cores=2, x_partition=2)
if __name__ == "__main__":
    main()