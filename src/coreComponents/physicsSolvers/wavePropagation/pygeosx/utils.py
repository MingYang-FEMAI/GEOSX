import json as j
import argparse


def parse_args():
    """Get parameters from python command and put them into args variable

    Return
    ------
    args : list
        Variable containing command parameters
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", help = "input xml file", required = True)
    parser.add_argument("--geosx", help = "GEOSX folder location", required = True)
    args = parser.parse_args()

    return args



def obj_to_dict(obj):
    if not  hasattr(obj,"__dict__"):
        return obj
    result = {}
    for key, val in obj.__dict__.items():
        if key.startswith("_"):
            continue
        element = []
        if isinstance(val, list):
            for item in val:
                element.append(todict(item))
        else:
            element = todict(val)
        result[key] = element
    return result


def ricker(maxT, dt, f0):
    """ Source function

    Parameters
    ----------
    maxT : float
        The max time for simulation

    dt : float
        The time step for simulation

    f0 : float
        Intensity

    Return
    ------
    fi :
        np array containing source value at all timestep
    """

    T0 = 1.0/f0;
    fi = np.zeros(int(maxT/dt))

    for t in range(int(maxT/dt)):
        t0 = dt*t
        if t0 <= -0.9*T0 or t0 >= 2.9*T0:
            fi[t] = 0.0;
        else:
            tmp      = f0*t0-1.0
            f0tm1_2  = 2*(tmp*np.pi)*(tmp*np.pi)
            gaussian = m.exp( -f0tm1_2 )
            fi[t]    = -(t0-1)*gaussian

    return fi


 def obj_to_json(obj,
                 filename=None):
     dic = obj_to_dict(obj)
     json = j.dumps(dic, indent=4)
     with open(filename + '.json', 'w') as jfile:
         print(json, end="", file=jfile)
         jfile.close()

     return filename + '.json'


def json_to_obj(jfile):
    obj = j.loads(studentJsonData, object_hook=lambda d: Namespace(**d))
    return obj