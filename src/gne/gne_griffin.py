import numpy as np
from astropy.cosmology import FlatLambdaCDM
from typing import Optional, Dict, Literal

ModelKind = Literal["shark", "galform"]

# Physical constants
C = 2.99792458e10        # cm/s
C2=C**2
MSUN = 1.98847e33        # g
YR = 3.15576e7           # s
G_PER_S = MSUN / YR      # (Msun/yr) -> g/s

def r_lso(spin):
        """ISCO radius in units of RG for Kerr BH (Bardeen et al. 1972) used in eqs. (3)-(6)."""
        spin = np.clip(np.asarray(spin, dtype=float), -0.999, 0.999)  # avoid singular edges
        abs_spin = np.abs(spin)
        Z1 = 1 + (1 - abs_spin**2)**(1/3) * ((1 + abs_spin)**(1/3) + (1 - abs_spin)**(1/3))
        Z2 = np.sqrt(3*abs_spin**2 + Z1**2)
        # prograde if spin>0, retrograde if spin<0
        r_lso = 3 + Z2 - np.sign(spin) * np.sqrt((3 - Z1) * (3 + Z1 + 2*Z2))
        return r_lso

def eps_thin(spin):
        """Thin-disc radiative efficiency ε_TD(a) from eq. (3)."""
        r = r_lso(spin)
        return 1 - np.sqrt(1 - 2/(3*r))

def L_edd(M_bh):
        """Eddington luminosity (erg/s)."""
        return 1.26e46 * (M_bh / 1e8)

def Mdot_edd(M_bh):
        """Eddington mass accretion rate (g/s) using the paper's fixed 0.1 efficiency convention."""
        return L_edd(M_bh) / (0.1 * C2)

def beta_from_alpha(alpha_ADAF=0.1):
        """Magnetic-to-total pressure parameterization used in ADAF formulae."""
        return 1 - alpha_ADAF/0.55  # Hawley et al. (1995) relation adopted in the paper

def luminosity_from_mdot(M_bh_msun,
        mdot_msun_per_yr, 
        spin=None,                     
        alpha_ADAF=0.1,
        delta_ADAF=0.2,
        eta_Edd=4.0,
        mcrit_ADAF=0.01
    ):
        """
        Compute L_bol (erg/s) for given BH mass (Msun), the accretion rate (Msun/yr), 
        and the dimensionless spin parameter (between -1 and 1).
        following eqs. (29)-(31) and (32). Vectorized.
        """
        if np.shape(M_bh_msun) != np.shape(mdot_msun_per_yr):
              raise ValueError(f"M_bh shape {np.shape(M_bh_msun)} != Mdot shape {np.shape(mdot_msun_per_yr)}")

        
        M_bh = np.asarray(M_bh_msun, dtype=float)
        Mdot = np.asarray(mdot_msun_per_yr, dtype=float) * G_PER_S  # g/s
        
        # Dimensionless mdot ≡ Mdot / Mdot_Edd with Mdot_Edd defined using 0.1 efficiency
        mdot_dim = Mdot / Mdot_edd(M_bh)

        # ε_TD from spin if provided; otherwise use 0.1 nominal
        if spin is None:
            epsTD =0.1
            r_lso_use = 6.0
        else:
           #clip the spin values to stay between -1 and 1 (physical values)
           spin = np.clip(np.asarray(spin, dtype=float), -0.999, 0.999)
           r_lso_use = r_lso(spin) 
           epsTD = eps_thin(spin)
            
        #I need to mask per case, but if no spin r_Iso and eps are a scalar the mask breaks    
        epsTD     = np.broadcast_to(epsTD,     np.shape(mdot_dim))
        r_lso_use = np.broadcast_to(r_lso_use, np.shape(mdot_dim))
        
        beta = beta_from_alpha(alpha_ADAF)
        # eq. (32): boundary inside ADAF so L is continuous
        mcrit_visc = 1e-3 * (delta_ADAF/5e-4) * ((1 - beta)/beta) * (alpha_ADAF**2)
        #New criteria from Griffin 2020 correction
        mcrit_super = eta_Edd * (0.1 / epsTD)           
        # Initialize Lbol
        Lbol = np.zeros_like(Mdot, dtype=float)

        # CASE 1: low-Ṁ ADAF (mdot < mcrit_visc) -> eq. (29)
        mask1 = mdot_dim < mcrit_visc
        if np.any(mask1):
            Lbol[mask1] = (
                2e-4 * epsTD[mask1] * Mdot[mask1] * C2
                * (delta_ADAF/5e-4)
                * ((1 - beta)/0.5)
                * (6.0 / r_lso_use[mask1])
            )

        # CASE 2: high-Ṁ ADAF (mcrit_visc ≤ mdot < mcrit_ADAF) -> eq. (30)
        mask2 = (mdot_dim >= mcrit_visc) & (mdot_dim < mcrit_ADAF)
        if np.any(mask2):
            Lbol[mask2] = (
                0.2 * epsTD[mask2] * Mdot[mask2] * C2
                * (mdot_dim[mask2] / (alpha_ADAF**2))
                * (beta/0.5)
                * (6.0 / r_lso_use[mask2])
            )

        # CASE 3: thin disc (mcrit_ADAF ≤ mdot ≤ mcrit_super) -> eq. (28)
        #mask3 = (mdot_dim >= mcrit_ADAF) & (mdot_dim <= eta_Edd) <--- old limit before 2020 correction
        mask3 = (mdot_dim >= mcrit_ADAF) & (mdot_dim <= mcrit_super)
        if np.any(mask3):
            Lbol[mask3] = epsTD[mask3] * Mdot[mask3] * C2

        # CASE 4: super-Eddington (mdot > mcrit_super) -> Griffin+2020 eq. (2)
        #mask4 = mdot_dim > eta_Edd<--- old limit before 2020 correction
        mask4 = mdot_dim > mcrit_super
        if np.any(mask4):
            arg = mdot_dim[mask4]/mcrit_super[mask4] #>1 by construction so log should be ok
            Lbol[mask4] = eta_Edd * (1 + np.log(arg)) * L_edd(M_bh[mask4])

        return Lbol
        

class BoolLuminosityFunction:
    """
    Compute luminosity functions (LFs) and starburst weights using
    a fixed snapshot's galaxy data + cosmology.
    """

    def __init__(self,
                 data: Dict[str, np.ndarray],
                 cosmology: Dict[str, float],
                 Lbox_Mpch: float,
                 n_batches: int,
                 tau_fold: float,
                 z_previous: float,
                 tot_batches: Optional[int] = None
                ):
        """
        data:
          For SHARK must include 'm_bulge','r_bulge' (r_bulge in comoving Mpc).
          For GALFORM must include 't_Q' (years).
        cosmology: dict with 'h','Omega_m','z'
        """
        self.data = data
        self.h = float(cosmology["h"])
        self.Om = float(cosmology["Omega_m"])
        self.z = float(cosmology["z"])
        self.Lbox = float(Lbox_Mpch)              # Mpc/h (comoving)
        self.n_batches = int(n_batches)
        self.tot_batches = int(tot_batches) if tot_batches is not None else int(n_batches)
        self.previous_z=float(z_previous)
        self.tau_fold=tau_fold
        # Volume in (Mpc)^3 to match φ units (h/Mpc)^3 dex^-1
        self.volume = (self.Lbox/self.h) ** 3
        self._cosmo_model = FlatLambdaCDM(H0=100.0 * self.h, Om0=self.Om)

    def compute_tQ(self) -> np.ndarray:
        """
        Return t_Q (years) for each galaxy.
       
        t_Q = fq * t_dyn, with
            r_phys_kpc = r_bulge_cMpc * 1000 / (1+z),
            M_b = m_bulge, Note that even if r is the half mass radious shark consideres 
            this m_bulge in the dynamic time computation so no 0.5 missing. 
            t_dyn = sqrt(r^3 / (G M_b)) [kpc/(km/s)] -> years.
        """
        # --- SHARK computation ---
        M_bulge = np.asarray(self.data["mstars_bulge"], dtype=float)       # Msun
        r_bulge_cMpc = np.asarray(self.data["rgas_bulge"], dtype=float)  # comoving Mpc

        # constants
        G = 4.30091e-6          # (kpc/Msun) (km/s)^2
        KMPCS = 3.085677581e16  # km per kpc
        SECYR = 3.15576e7       # s per yr

        t_bulge_yr = np.zeros_like(M_bulge, dtype=float)
        mask = (
            np.isfinite(M_bulge) & (M_bulge > 0) &
            np.isfinite(r_bulge_cMpc) & (r_bulge_cMpc > 0)
        )
        if np.any(mask):
            r_phys_kpc = (r_bulge_cMpc[mask] * 1000.0) / (1.0 + self.z)
            t_dyn_kpc_per_kms = np.sqrt(r_phys_kpc**3 / (G * M_bulge[mask]))
            t_dyn_sec = t_dyn_kpc_per_kms * KMPCS
            t_bulge_yr[mask] = t_dyn_sec / SECYR

        t_Q_yr = self.tau_fold * t_bulge_yr
        return t_Q_yr

    def compute_sb_weights(self) -> np.ndarray:
        """
        Starburst duty weights:
            w_sb = min(t_Q / Δt, 1)
        where Δt is the time between previous_z and z.
        """
        t_Q_yr = self.compute_tQ()
        dt_window_yr = (self._cosmo_model.age(self.z).value-self._cosmo_model.age(self.previous_z).value) * 1e9
        return np.minimum(t_Q_yr / dt_window_yr, 1.0)

    def blf(self,
           L_bol: np.ndarray,
           dp: float = 0.3,
           edges: Optional[np.ndarray] = None,
           compute_weights: bool = True):
        """
        Differential LF:
            φ(log10 L) = dN / (dV d log10 L)  [ (h/Mpc)^3 dex^-1 ]
        Normalization: counts * (tot_batches / n_batches) / (volume * dp)
        """
        L_bol = np.asarray(L_bol, dtype=float)
        mask = np.isfinite(L_bol) & (L_bol > 0.0)
        if not np.any(mask):
            raise ValueError("No valid positive luminosities provided.")

        logL = np.log10(L_bol[mask])
        if edges is None:
            edges = np.arange(42.0, 49.0 + dp * 0.5, dp)
        centers = 0.5 * (edges[1:] + edges[:-1])

        if compute_weights:
            w = self.compute_sb_weights()
            w = w[mask]
            counts, _ = np.histogram(logL, bins=edges, weights=w)

        else:
            counts, _ = np.histogram(logL, bins=edges)
        phi = counts * (self.tot_batches / self.n_batches) / (self.volume * dp)
        return centers, phi

    def instantaneous_luminosity(self, L_sb: np.ndarray):
        # 1) Bernoulli draw: SB is “on” with probability p_on
        weights = self.compute_sb_weights()
        rng = np.random.default_rng(42)              
        on_sb = rng.random(weights.shape[0]) < weights   
        # 2) Instantaneous SB luminosity (0 if off)
        L_inst_sb = np.where(on_sb, L_sb, 0.0)   # erg/s
        return L_inst_sb

def compute_Lbol_griffin(data,h0,omega0,redshift,Lbox,params,redshift_previous,N_batches=1,tau_fold=None,tot_batches=1):
    
    if tau_fold is None:
        tau_fold = 1

    data_dict = {params[i_param].split('/')[-1]: data[i_param] for i_param in range(len(params))}

    cosmo = {'h': h0, 'Omega_m': omega0, 'z': redshift}
    print('Computing luminosity for the sb component')
    blf = BoolLuminosityFunction(data=data_dict,cosmology=cosmo,Lbox_Mpch=Lbox,
                                 n_batches=N_batches,tot_batches=tot_batches,
                                 tau_fold=tau_fold,z_previous=redshift_previous
                                 )
    bol_luminosity_sb_insnapshot=luminosity_from_mdot(data_dict["m_bh"],data_dict["bh_accretion_rate_sb"])
    bol_luminosity_sb=blf.instantaneous_luminosity(bol_luminosity_sb_insnapshot)
    print('Computing luminosity for the hh component')
    bol_luminosity_hh=luminosity_from_mdot(data_dict["m_bh"],data_dict["bh_accretion_rate_hh"])
    return bol_luminosity_hh+bol_luminosity_sb
    