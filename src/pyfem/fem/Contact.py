
import numpy as np
from numpy import append, repeat


# ===============================================================================
#
# ===============================================================================


class Contact:

    def __init__(self, props):

        self.flag = False
        self.dispDofs = ["u", "v"]
        self.centre = [1., 1.]
        self.direction = [0.0, 0.0]
        self.radius = 10.
        self.penalty = 1.e6

        if hasattr(props, 'contact'):
            if hasattr(props.contact, 'type'):
                self.type = props.contact.type
                self.flag = True
                self.centre = np.array(props.contact.centre)
                self.radius = props.contact.radius
                self.penalty = props.contact.penalty
                self.direction = np.array(props.contact.direction)

    # -------------------------------------------------------------------------------
    #  checkContact   (with flag)
    # -------------------------------------------------------------------------------

    def checkContact(self, row, val, col, B, globdat):

        if not self.flag:
            return

        centre = self.centre + globdat.lam * self.direction

        for node_id in list(globdat.nodes.keys()):
            crd = globdat.nodes.get_node_coords(node_id)

            idofs = globdat.dofs.get_dof_ids_by_types([node_id], self.dispDofs)

            crd += globdat.state[idofs]

            ds = crd - centre

            dsnorm = np.linalg.norm(ds)
            overlap = self.radius - dsnorm

            if overlap > 0:
                normal = ds / dsnorm

                B[idofs] += -self.penalty * overlap * normal

                mat = self.penalty * np.outer(normal, normal)

                row = append(row, repeat(idofs, len(idofs)))

                for i in range(len(idofs)):
                    col = append(col, idofs)

                val = append(val, mat.reshape(len(idofs) * len(idofs)))

        return row, val, col
