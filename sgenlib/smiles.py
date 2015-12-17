# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Routines to manipulate with SMILES
"""

import urllib

class RetrieveSmilesException(Exception):
    pass

def _get_cactus(smiles):
    """
    Convert SMILES using Cactus service

    Parameters
    ----------
    smiles : string
      the SMILES string

    Returns
    -------
    string
      the path to the downloaded SDF file

    Raises
    ------
    RetrieveSmilesException
      if the retrieval fails
    """
    site = "cactus.nci.nih.gov"
    try:
        page = urllib.urlopen("http://%s/cgi-bin/translate.tcl?smiles=%s&format=sdf&astyle=kekule&dim=3D&file="%(site, smiles))
    except:
        raise RetrieveSmilesException("Could not connect with server")

    for line in page:
        if "Click here" in line  and 'a href="' in line :
            dummy1, url, dummy2 = line.split('"')
            try:
                path, header = urllib.urlretrieve("http://%s%s"%(site,url))
            except:
                raise RetrieveSmilesException("Could not retrieve file")
            return path


def _get_indiana(smiles) :
    """
    Convert SMILES using Indiana service

    Parameters
    ----------
    smiles : string
      the SMILES string

    Returns
    -------
    string
       the path to the downloaded SDF file

    Raises
    ------
    RetrieveSmilesException
      if the retrieval fails
    """
    try :
        path, header = urllib.urlretrieve("http://cheminfov.informatics.indiana.edu/rest/thread/d3.py/SMILES/%s" % smiles)
    except :
        raise RetrieveSmilesException("Could not connect with server")
    else:
        with open(path,"r") as f:
            if f.readline().startswith("<!DOCTYPE HTML"):
                raise RetrieveSmilesException("Server Error")
        return path

def _sdf2xyz(sdffile,xyzfile) :
    """
    Read a SDF file and convert it to xyz-format

    Parameters
    ----------
    sdffile : string
      the name of the SDF file
    xyzfile : string
      the name of the XYZ-file
    """
    with open(sdffile,"r") as f :
        lines = f.readlines()

    if len(lines) == 0 : return

    with open(xyzfile,"w") as f :
        natom = int(lines[3].strip().split()[0])
        f.write("%d\n0\n"%natom)
        for line in lines[4:4+natom] :
            cols = line.strip().split()
            f.write("%s %s %s %s\n"%(cols[3],cols[0],cols[1],cols[2]))


def convert2xyz(smiles,out,web=None,verbose=True):
    """
    Convert a SMILES string to an xyz file using either of two web services

    Parameters
    ----------
    smiles : string
      the SMILES strings
    out : string
      the name of the xyz file
    web : string, optional
      can be either cactus, indiana, both or none
      determines the web service to use
    verbose : boolean, optional
      if the routine should write out verbose information
    """

    # Try to convert the SMILES string using one of two web services
    service = None
    if web == "cactus" :
        try:
            path = _get_cactus(smiles.strip())
        except:
            if verbose :
                print "Unable to convert SMILES."
            return
        else:
            service = "cactus"
    elif web == "indiana" :
        try:
            path = _get_indiana(smiles.strip())
        except:
            if verbose :
                print "Unable to convert SMILES."
            return
        else:
            serivce = "indiana"
    else :
        try:
            path = _get_indiana(smiles.strip())
        except:
            try:
                path = _get_cactus(smiles.strip())
            except:
                if verbose :
                    print "Unable to convert SMILES."
                return
            else:
                service = "cactus"
        else:
            service = "indiana"


    text = {"indiana":"the smi23d web service provided by the Chemical Informatics and Cyberinfrastructure Collaboratory at Indiana University","cactus":" the SMILES translator provided by the National Cancer Institute CADD group"}
    link = {"indiana":"http://www.chembiogrid.info/projects/proj_ws_all.html","cactus":"http://cactus.nci.nih.gov/translate/"}
    _sdf2xyz(path,out)
    if verbose :
        print "SMILES converted sucessfully using %s"%text[service]
        print "Web page: %s"%link[service]

    return out
