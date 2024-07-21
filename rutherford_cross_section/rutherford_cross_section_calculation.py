import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use('seaborn-poster')

# Set the constant in function to the proyectile and target
def sigma_laboratory(ec=1000, zp=1, zt=6, ap=1, at=12, theta_deg=170):
    # SIMNRA formnula for calculation of Rutherford differential cross section in Lab coordinates
    # simga [mb/sr], Energy [keV]
    theta_lab = np.radians(theta_deg)
    const = 5.1837436 * 10 ** 6
    term1 = (zp * zt / ec) ** 2
    root_term = np.sqrt( ( at ** 2 ) - ( ap * ap * np.sin(theta_lab) * np.sin(theta_lab) ) )
    sigma = const * term1 * ( (root_term + at * np.cos(theta_lab)) ** 2 ) / ( at * np.sin(theta_lab) ** 4 * root_term )
    return round(sigma,3)
    #return root_term, np.sin(np.radians(170))
    

def sigma_center_of_mass(ec=1, zp=1, zt=13, theta_deg=170):
    # Rutherford cross section in Center of Mass coordinates frame
    # simga [mb/sr], Energy [MeV]
    theta_cm_rad = np.radians(theta_deg)
    sigma = 1.296 * (zp * zt / ec) * (1 / np.sin(theta_cm_rad / 2) ** 4)
    return round(sigma,3) # falta revisar que de resultados correctos


def convert_sigma_lab_to_cm(sigma_lab, ap=1, at=27, theta_deg=170):
    # Convertion from diff cross section of Laboratory to Center of Mass coordinates
    theta_lab = np.radians(theta_deg)
    lan = ap / at
    theta_cm = theta_lab + np.arcsin(lan * np.sin(theta_lab))
    jacobian = ( 1 + lan * np.cos(theta_lab) ) / ( (1 + lan * lan + 2 * lan * np.cos(theta_lab)) ** (3/2) )
    sigma_cm = sigma_lab * jacobian
    return round(sigma_cm,3)

def open_data(name='data_12c_p.dat'):
    filedir = os.getcwd()
    denergy = np.loadtxt(filedir+"/"+name, usecols=(0,), skiprows=1)
    dsigma = np.loadtxt(filedir+"/"+name, usecols=(1,), skiprows=1)
    return denergy, dsigma

def plot(exp_energy, exp_sigma, calc_energy, calc_sigma):
    filedir = os.getcwd()
    
    fig, ax = plt.subplots(1,1,figsize=(10, 6)) 

    ax.plot(exp_energy, exp_sigma, '.',color='blue', lw=1.5, label=r'sección eficaz experimental')
    ax.plot(calc_energy, calc_sigma, '--',color='red', lw=1.5, label=r'sección eficaz de Rutherford')


    ax.grid(True, color='gray', linestyle='--', linewidth=0.5, which='both', alpha=0.5)
    #ax.set(xlim=(500, 2400))
    #ax.set(ylim=(0, 1200))
    ax.set_xlabel(r'Energía [keV]')
    ax.set_ylabel(r'Sección Eficaz de Rutherford [mb/sr]')
    ax.set_title('$^{12}$C(p,p0)')
    ax.legend(loc='best', prop={'size': 10})

    plt.tight_layout()
    plt.savefig(filedir+"/seccion_carbono.png", bbox_inches="tight", dpi=300)
    plt.show()
   


if __name__ == "__main__":
    # Energy
    lst_energy = np.arange(400, 2410, 10) # energy in keV
    sigma = list(map(sigma_laboratory, lst_energy))
    sigma_calc = np.array(sigma)

    de, dsig = open_data()

    plot(de, dsig, lst_energy, sigma_calc)

    #print(lst_energy,sigma_calc)
    #print(lst_energy)

    #print(type(sigma_calc))    
    #print(sigma_laboratory())
    #print(convert_sigma_lab_to_cm(sigma_laboratory()))
    