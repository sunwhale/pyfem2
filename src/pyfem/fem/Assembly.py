from numpy import zeros, ones, append, repeat, array
from scipy.sparse import coo_matrix

from pyfem.utils.data_structures import elementData


def assembleArray(props, globdat, rank, action):
    # Initialize the global array A with rank 2

    B = zeros(globdat.dofs.get_number_of_dofs() * ones(1, dtype=int))
    cc = 0.0

    val = array([], dtype=float)
    row = array([], dtype=int)
    col = array([], dtype=int)

    nDof = globdat.dofs.get_number_of_dofs()

    if action != 'commit':
        globdat.resetNodalOutput()

    # Loop over the element groups
    for elementGroup in globdat.elements.iter_group_names():

        # Get the properties corresponding to the elementGroup
        el_props = getattr(props, elementGroup)

        # Loop over the elements in the elementGroup
        for iElm, element in enumerate(globdat.elements.iter_element_group(elementGroup)):

            # Get the element nodes
            el_nodes = element.getNodes()

            # Get the element coordinates
            el_coords = globdat.nodes.get_node_coords(el_nodes)

            # Get the element degrees of freedom
            el_dofs = globdat.dofs.get_dof_ids_by_types(el_nodes, element.dof_types)

            # Get the element state
            el_a = globdat.state[el_dofs]
            el_Da = globdat.dstate[el_dofs]

            # Create the an element state to pass through to the element
            # el_state = Properties( { 'state' : el_a, 'dstate' : el_Da } )
            elemdat = elementData(el_a, el_Da)

            elemdat.coords = el_coords
            elemdat.nodes = el_nodes
            elemdat.props = el_props
            elemdat.iElm = iElm

            element.globdat = globdat

            if hasattr(element, "matProps"):
                elemdat.matprops = element.matProps

            if hasattr(element, "mat"):
                element.mat.reset()

            # Get the element contribution by calling the specified action
            if hasattr(element, action):
                getattr(element, action)(elemdat)

            # for label in elemdat.outlabel:
            #  element.appendNodalOutput( label , globdat , elemdat.outdata )

            # Assemble in the global array
            if rank == 1:
                B[el_dofs] += elemdat.fint
                cc += elemdat.diss
            elif rank == 2 and action == "getTangentStiffness":

                row = append(row, repeat(el_dofs, len(el_dofs)))

                for i in range(len(el_dofs)):
                    col = append(col, el_dofs)

                val = append(val, elemdat.stiff.reshape(len(el_dofs) * len(el_dofs)))

                B[el_dofs] += elemdat.fint
            elif rank == 2 and action == "getMassMatrix":

                row = append(row, repeat(el_dofs, len(el_dofs)))

                for i in range(len(el_dofs)):
                    col = append(col, el_dofs)

                val = append(val, elemdat.mass.reshape(len(el_dofs) * len(el_dofs)))

                B[el_dofs] += elemdat.lumped
    #    else:
    #      raise NotImplementedError('assemleArray is only implemented for vectors and matrices.')

    if rank == 1:
        return B, cc
    elif rank == 2:

        '''if globdat.contact.flag:
          row , val , col = globdat.contact.checkContact( row , val , col , B , globdat )      '''

        return coo_matrix((val, (row, col)), shape=(nDof, nDof)), B


def assembleInternalForce(props, globdat):
    fint = assembleArray(props, globdat, rank=1, action='getInternalForce')
    return fint[0]


def assembleExternalForce(props, globdat):
    fext = assembleArray(props, globdat, rank=1, action='getExternalForce')

    return fext[0] + globdat.fhat * globdat.solver_status.lam


def assembleDissipation(props, globdat):
    return assembleArray(props, globdat, rank=1, action='getDissipation')


def assembleTangentStiffness(props, globdat):
    return assembleArray(props, globdat, rank=2, action='getTangentStiffness')


def assembleMassMatrix(props, globdat):
    return assembleArray(props, globdat, rank=2, action='getMassMatrix')


def commit(props, globdat):
    return assembleArray(props, globdat, rank=0, action='commit')


def getAllConstraints(props, globdat):
    # Loop over the element groups
    for elementGroup in globdat.elements.iter_group_names():

        # Get the properties corresponding to the elementGroup
        el_props = getattr(props, elementGroup)

        # Loop over the elements in the elementGroup
        for element in globdat.elements.iter_element_group(elementGroup):
            # Get the element nodes
            el_nodes = element.getNodes()

            elemdat.nodes = el_nodes
            elemdat.props = el_props

            # Get the element contribution by calling the specified action
            getattr(element, 'getConstraints', None)(elemdat)
