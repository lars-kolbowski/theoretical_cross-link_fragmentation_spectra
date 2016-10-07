from pyteomics import mass
import re
import numpy as np
import pandas as pd
import os

class MS2_spectrum():
    """
    Class container for MS2 spectra.
    We need the following input parameters:
    title, RT, pepmass, pepint, charge, peaks, peakcharge=[]

    Parameters:
    -----------------------------------------
    title: str,
            title of the spectrum
    RT: float,
        retention time of the precursor
    pepmass: float,
              mass of the precursor
    charge: int,
             charge of the precursor
    peaks: [(float, float)],
           mass intensity
    peakcharge: arr,
                charge array for the peaks

    """
    def __init__(self, title, RT, pepmass, pepint, charge, peaks, peakcharge=[]):
        self.title = title
        self.RT = RT
        self.pepmass = pepmass
        self.pepint = pepint
        self.charge = charge
        self.peaks = peaks
        self.peakcharge = peakcharge

    def getPrecursorMass(self):
        """
        Returns the precursor mass
        """
        return(self.pepmass)

    def getPrecursorIntensity(self):
        """
        Returns the precursor mass
        """
        return(self.pepint)

    def getRT(self):
        """
        Returns the precursor mass
        """
        return(self.RT)

    def getTitle(self):
        """
        Returns the precursor mass
        """
        return(self.title)

    def getPeaks(self):
        """
        Returns the precursor mass
        """
        return(self.peaks)

    def getMasses(self):
        """
        TODO:
        Returns the precursor mass
        """
        return(self.peaks[:,0])

    def getIntensities(self):
        """
        TODO:
        Returns the precursor mass
        """
        return(self.peaks[:,1])

    def getUnchargedMass(self):
        """
        Computs the uncharged mass of a fragment:
        uncharged_mass = (mz * z ) - z
        """
        return( (self.pepmass * self.charge) -  self.charge)

    def printf(self):
        print ("Title, RT, PEPMASS, PEPINT, CHARGE")
        print (self.title, self.RT, self.pepmass, self.pepint, self.charge)

    def to_mgf(self):
        # need dummy values in case no peak charges are in the data
        if len(self.peakcharge) == 0:
            self.peakcharge = [""]*self.peaks.shape[0]
        mgf_str="""
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (self.title, self.RT, self.pepmass, self.pepint, self.charge, "\r\n".join(["%s %s %s" % (i[0], i[1], j, ) for i,j in zip(self.peaks, self.peakcharge)]))
        return(mgf_str)

    def to_mgf_Xi(self):
        # need dummy values in case no peak charges are in the data
        if len(self.peakcharge) == 0:
            self.peakcharge = [""]*self.peaks.shape[0]
        mgf_str="""
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (self.title, self.RT, self.pepmass, self.pepint, self.charge, "\r\n".join(["%s %s" % (i[0], i[1]) for i in self.peaks]))
        return(mgf_str)


def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxcharge`.
    """
    for i in xrange(1, len(peptide)-1):
        for ion_type in types:
            for charge in xrange(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)


def extract_mods(peptide):
    #        """
    #        Function that extracts the modification info from the peptide seq.
    #        returns the modification info (name, mass, position) in a list
    #
    #        Parameters:
    #        -----------------------------
    #        peptide: str,
    #                    peptide sequence
    #
    #        Return: modification info with masses and position
    #        """

    mods_mass_list = [
        {'name': 'bs3', 'mass': 27.984},
        {'name': 'bs3loop', 'mass': 138.06807},
        {'name': 'bs3nh2', 'mass': 155.094619105},
        {'name': 'bs3oh', 'mass': 156.0786347},
        {'name': 'cm', 'mass': 57.021464},
        {'name': 'ox', 'mass': 15.994915},
    ]

    mods = []
    offset = 0
    for mod in [m for m in re.finditer('[a-z0-9]+', peptide)]:

        # get modification mass from mods_mass_list
        found = False
        for m in mods_mass_list:
            if mod.group() == m['name']:
                mod_mass = m['mass']
                found = True
        if not found:
            print mod.group()
        # get position of modification
        pos = mod.start() + offset

        mods.append({'name': mod.group(), 'pos': pos, 'mass': mod_mass})
        offset += mod.start() - mod.end()

    return mods


def filter_out_charges(fragment, pep1, pep2=""):
    peptide = fragment[0]
    theoretical_charge = 0
    if "+P" in fragment[0]:
        theoretical_charge += 1
        peptide = peptide[:-2]
        # if it's a CL peptide filter out charge state 1
        if fragment[2] == 1:
            return False
        # intensity 100 means its a pep1 fragment - 50 pep2 fragment
        if fragment[3] == 100:
            peptide += pep2
        else:
            peptide += pep1
    theoretical_charge += sum([1 for l in peptide if l in "RNQHK"])
    if len(peptide) > theoretical_charge:
        theoretical_charge += 1

    if theoretical_charge < fragment[2]:
        return False
    else:
        return True


def CL_fragments(peptide, modifications, mod_mass, link=-1, types=('b', 'y'), maxcharge=1, annotation="fragments"):
    #         """
    #        The function generates all possible m/z for fragments of types
    #        `types` and of charges from 1 to `maxcharge`.
    #
    #        Parameters:
    #        ------------------------------
    #        peptide: str,
    #                peptide sequence
    #        mod_mass: float,
    #                    modification mass (e.g. the mass of the cross-linked peptide)
    #        link: int,
    #                the index where the cross-link or modification is located
    #        """
    #print "XLINK at:", peptide[link]
    #yindex = len(peptide) - link

    for i in xrange(1, len(peptide)):
        for ion_type in types:
            for charge in xrange(1, maxcharge+1):
                #abc ions
                if ion_type[0] in 'abc':
                    if annotation == "ions":
                        res_name = ion_type[0] + str(i)
                    else:
                        res_name = peptide[:i]
                    res_mass = mass.fast_mass(peptide[:i], ion_type=ion_type, charge=charge)
                    offset = 0
                    for mod in modifications:
                        if i >= mod['pos']:
                            res_mass += mod['mass']/charge
                            if not annotation == "ions":
                                res_name = res_name[:mod['pos']+offset] + mod['name'] + res_name[mod['pos']+offset:]
                            offset += len(mod['name'])
                    if i >= link and not link == -1:
                        res_name += "+P"
                        res_mass += mod_mass/charge
                    yield [res_name, res_mass, charge]

                # xyz
                else:
                    if annotation == "ions":
                        res_name = ion_type[0] + str(len(peptide)-i)
                    else:
                        res_name = peptide[i:]
                    res_mass = mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge)
                    offset = 0
                    for mod in modifications:
                        if i < mod['pos']:
                            res_mass += mod['mass']/charge
                            if not annotation == "ions":
                                res_name = res_name[:mod['pos']-i+offset] + mod['name'] + res_name[mod['pos']-i+offset:]
                            offset += len(mod['name'])
                    if i < link and not link == -1:
                        res_name += "+P"
                        res_mass += mod_mass/charge
                    yield [res_name, res_mass, charge]


def get_theoretical_spectrum(peptide1, peptide2="", cl_mass=138.06807961, link1=0, link2=0, types=('b', 'y'), maxcharge=1, annotation="fragments"):
    #        """
    #        Function that generates theoretical spectra. It is capable to generate
    #        cross-linked spectra as well.
    #
    #        Parameters:
    #        -----------------------------
    #        peptide1: str,
    #                    peptide sequence
    #        peptide2: str,
    #                    peptide sequence 2
    #        cl: bool,
    #           cross-linked (True) or not (False)
    #         link1, link2: int,
    #                 position of the respective cross-linking site for peptide 1/2
    #
    #        """

    mods1 = extract_mods(peptide1)
    # delete mods from AA seq
    peptide1 = re.sub("[a-z0-9]+", "", peptide1)

    if peptide2 == "":
        pep1mass = mass.fast_mass(peptide1) + sum([m['mass'] for m in mods1])
        masses = [i for i in CL_fragments(peptide1, mods1, pep1mass, types=types, maxcharge=maxcharge, annotation=annotation)]
        for m in masses:
            m.append(100)
        the_fragments = masses
        #the_fragments = sorted([i for i in fragments(peptide1)])
        filtered_fragments = the_fragments
    else:
        mods2 = extract_mods(peptide2)
        # delete mods from AA seq
        peptide2 = re.sub("[a-z0-9]+", "", peptide2)

        pep1mass = mass.fast_mass(peptide1) + sum([m['mass'] for m in mods1])
        pep2mass = mass.fast_mass(peptide2) + sum([m['mass'] for m in mods2])
        masses_1 = [i for i in CL_fragments(peptide1, mods1, pep2mass+cl_mass, link1, types=types, maxcharge=maxcharge)]
        for m in masses_1:
            m.append(100)
        masses_2 = [i for i in CL_fragments(peptide2, mods2,  pep1mass+cl_mass, link2, types=types, maxcharge=maxcharge)]
        for m in masses_2:
            m.append(50)
        the_fragments = masses_1 + masses_2
        precursor_mass = mass.fast_mass(peptide1) + mass.fast_mass(peptide2) + cl_mass

        #find maxcharge for fragment depending on AAseq
        filtered_fragments = [f for f in the_fragments if filter_out_charges(f, peptide1, peptide2)]


    return filtered_fragments


#get_theoretical_spectrum("AACcmLLPKLDELRDEGKASSAKQR", "LKCcmASLQK", link1=21, link2=2, types=('b', 'c', 'x', 'y'), maxcharge=3)
#get_theoretical_spectrum("AEFAEVSKLVTDLTK", "AFKAWAVAR", link1=8, link2=3)


#####################################
def save_theoretical_spectra_as_mgf(pepInfosCSV, protein):
    df = pd.read_csv(pepInfosCSV)
    #df is the dataframe with the input peptides that should be simulated
    df.columns.values.tolist()

    method_ions = [
        {'method': 'ETD', 'ions': ('c', 'z')},
        {'method': 'ETciD', 'ions': ('b', 'c', 'y', 'z')},
        {'method': 'EThcD', 'ions': ('b', 'c', 'y', 'z')},
        {'method': 'CID', 'ions': ('b', 'y')},
        {'method': 'HCD', 'ions': ('b', 'y')},
    ]

    for method in method_ions:
        counter = 1
        scannumber = []
        spectra = []
        for row_i in df.iterrows():
            p1 = row_i[1].pep1.strip()
            p2 = row_i[1].pep2.strip()
            site1 = row_i[1].linkpos1   # + 1
            site2 = row_i[1].linkpos2   # + 1
            precursorMZ = row_i[1]['m/z']
            #convert peptide sequence to correct format
            #p1_long = p1.replace("cm", "(Carbamidomethyl)").replace("ox", "(Oxidation)")
            #p2_long = p2.replace("cm", "(Carbamidomethyl)").replace("ox", "(Oxidation)")
            spectrum = get_theoretical_spectrum(p1, p2, link1=site1, link2=site2, types=method['ions'], maxcharge=row_i[1].z)
            return spectrum
            mz = []
            inten = []
            charge = []
            fragment_names = []
            for fragment in spectrum:
                mz.append(fragment[1])
                inten.append(fragment[3])
                charge.append(fragment[2])
                fragment_names.append(fragment[0])
            #precursor = SpectrumGenerator.getPrecursorMZ(p1_long, p2_long, row_i[1].charge)
            peaks = np.array(zip(mz, inten))

            scannumber.append(str(counter).zfill(4))
            title = "RawFile:{}-{}_{}-{}_{} FinneganScanNumber: {} "\
                    .format(p1, site1, p2, site2, precursorMZ, str(counter).zfill(4))

            spectra.append(MS2_spectrum(title, 600, precursorMZ,
                                                         pepint=5000, charge=row_i[1].z,
                                                         peaks=peaks, peakcharge=fragment_names))
            counter += 1

        outfile = open("{}/{}simspectra.mgf".format("/home/lars/data/theoretical_spectra/"+protein+"/", method['method']), "w")
        for i in spectra:
            outfile.write(i.to_mgf().replace("\r\n", "\n"))
        outfile.close()

        xi_output_path = "/home/lars/data/theoretical_spectra/"+protein+"/Xi"
        if not os.path.exists(xi_output_path):
            os.makedirs(xi_output_path)
        outfile = open("{}/{}simspectra.mgf".format(xi_output_path, method['method']), "w")
        for i in spectra:
            outfile.write(i.to_mgf_Xi().replace("\r\n", "\n"))
        outfile.close()
