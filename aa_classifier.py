from topol.matrix_pro import proline
from topol.matrix_ala import alanine
from topol.matrix_gly import glycine

class AAClassifier():
    def __init__(self, aa_sequence):
        """
        Reads the primary sequence and returns the required angles.
        Can calculated the ramachandran bin id from given set of 
        angles. This class is a container for a specific sequence.
        """
        self.aa_sequence        = aa_sequence
        self.required_dihedrals = self.__generateRequiredDihedrals()
         
    def __getProlineId(self, omega, psi):
        """
        @return str(proline bin id)
        input angles range from [-180.0, 180.0[
        """
        omega = int(omega + 180)
        psi   = int(179 - psi)
        return str(proline[psi][omega])
    
    def __getAlanineId(self, phi, psi):
        """
        @return str(alanine bin id)
        """
        phi   = int(phi + 180)
        psi   = int(179 - psi)
        return str(alanine[psi][phi])
    
    def __getGlycineId(self, phi, psi):
        """
        @return str(alanine bin id)
        """
        phi   = int(phi + 180)
        psi   = int(179 - psi)
        return str(glycine[psi][phi])
    
    
    def getBinId(self, dihedrals):
        """
        @param dihedrals get from cpptraj
        returns the Bin Id
        """
        # Create a list of residue dicts with all the attributes
        # for easy access
        # dihedral_data = [{'res_id': 2, 'aa_type':'A', 'phi':value, 'psi':value<, 'omega':value>}, 
        #                  {}, ...]
        dihedral_data = []
        index = 0
        for residue in self.required_dihedrals:
            res_id       = residue['res_id']
            aa_type      = residue['aa_type']
            
            
            dihedral_data.append({'res_id'  : res_id,
                                  'aa_type' : aa_type})
            for angle in residue['required_angles']:
                dihedral_data[-1][angle] = dihedrals[index]
                index += 1
        
        # Generate the bin id tag now
        bin_id = ""
        for residue in dihedral_data:
            if residue['aa_type'] == 'P':
                bin_id += self.__getProlineId(residue['omega'], residue['psi'])
                continue
            elif residue['aa_type'] == 'G':
                bin_id += self.__getGlycineId(residue['phi'], residue['psi'])
                continue        
            else:
                bin_id += self.__getAlanineId(residue['phi'], residue['psi'])

        return bin_id

    def __generateRequiredDihedrals(self):
        """
        Should only be called once at startup
        @return [{'res_id'          : 2, 
                  'aa_type'         : 'A', 
                  'required_angles' : ['phi','psi']},
                 {...}, ...]
        """
        required_dihedrals = []
        for res_id, res_name in enumerate(self.aa_sequence):
            if   res_name in ['b' , 'e', '-']:
                continue
            elif res_name == 'P':
                required_dihedrals.append({'res_id'          : res_id + 1, 
                                           'aa_type'         : res_name, 
                                           'required_angles' : ['psi','omega']})
            else:
                required_dihedrals.append({'res_id'          : res_id + 1, 
                                           'aa_type'         : res_name, 
                                           'required_angles' : ['phi','psi']})
        return required_dihedrals
    
    def getRequiredDihedrals(self):
        return self.required_dihedrals
