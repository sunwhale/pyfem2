def dummy_print(something):
    r"""
    a print function which is the same as the original func

    :param str something: a str you want to print
    :return: no return

    .. math::
        {\varepsilon _{ij}} = \frac{1}{2}\left( {{u_{i,j}} + {u_{j,i}}} \right)
    .. math::
        \iiint\limits_V {{{\left( {\frac{{\partial A}}{{\partial {\varepsilon _{ij}}}}\delta {u_i}} \right)}_{,j}}{\text{d}}V} = \iint\limits_S {\frac{{\partial A}}{{\partial {\varepsilon _{ij}}}}\delta {u_i}{n_j}{\text{d}}S}
    
    这是一个简单的例子：

    Example
    -------
    >>> dummy_print("hello")
    hello
    """
    print(something)