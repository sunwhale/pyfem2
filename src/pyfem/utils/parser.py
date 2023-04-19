from pyfem.utils.data_structures import Properties


def has_value(db, val):
    """
    Addition of two numbers

    :param db: a
    :type db:  integer
    :param val: btted
    :type val:  integer
    :returns: new value
    :rtype:   integer
    """

    keys = list(db.keys())

    for key in keys:

        if type(db[key]) == dict:
            if has_value(db[key], val):
                return True
        else:
            if db[key] == val:
                return True
    return False


def clean_variable(a):
    """
    这个函数用于 "清理" 变量，它有一个输入参数 a，代表变量名或字符串值。

    函数的逻辑如下：

    如果 a 是字符串 "true"，则返回 True。
    如果 a 是字符串 "false"，则返回 False。
    否则，尝试将字符串转换为变量或表达式的值，如果转换成功，则返回这个值。
    如果转换失败或者 a 不是字符串，则返回 a 的原始值。
    
    :param a:
    :return:
    """
    if a == 'true':
        return True
    elif a == 'false':
        return False
    else:
        try:
            return eval(a)
        except:
            return a


def is_node_dof(node_dof):
    return type(node_dof) == str and '[' in node_dof


def decode_node_dof(node_dof, nodes):
    a = node_dof.split('[')
    dof_type = a[0]
    node_id = clean_variable(a[1].split(']')[0])

    if type(node_id) == str:
        return dof_type, nodes.groups[node_id]
    else:
        return dof_type, [node_id]


def get_type(a):
    return type(clean_variable(a))


def store_value(db, key, a):
    if type(a) == list:
        tmp = []
        for v in a:
            tmp.append(clean_variable(v))
        db.store(key, tmp)
    else:
        db.store(key, clean_variable(a))


def read_item(l1, db):
    if '.' in l1[0]:
        l2 = l1[0].split('.', 1)

        if l2[0] in db:
            if type(db[l2[0]]) == dict:
                child = db[l2[0]]
            else:
                child = Properties()
        else:
            child = Properties()

        l1[0] = l2[1]

        ln = read_item(l1, child)

        db[l2[0]] = child

        return ln

    else:
        l2 = l1[1].split(';', 1)

        if l2[0][0] == '[':
            l3 = l2[0][1:-1].split(',')
            store_value(db, l1[0], l3)
        else:
            store_value(db, l1[0], l2[0])

        return l2[1]


def read_block(ln, db):
    while True:

        if ln[0:7] == 'include':
            l1 = ln.split(';', 1)
            deep_file_name(l1[0][8:-1], db)
            ln = l1[1]
            continue

        l1 = ln.split('=', 1)

        if len(l1) == 1:
            return ln

        if l1[0][0:2] == '};':
            return ln[2:]

        if l1[0][0:2] == '//':
            ln = l1[1].split(';', 1)[1]
            continue

        # if l1[0][0:1] == '#':
        #  ln = l1[1].split(';',1)[1]
        #  continue

        if l1[1][0] == '{':
            child = Properties()
            ln = l1[1][1:]

            ln = read_block(ln, child)

            db.store(l1[0], child)

        else:
            ln = read_item(l1, db)


def file_parser(file_name):
    db = Properties()

    f = open(file_name)

    f2 = ''

    for line in f:
        if not line.startswith('#'):
            f2 = f2 + line

    # ln = open(file_name).read().replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '')
    ln = f2.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '')

    read_block(ln, db)

    return db


def deep_file_name(file_name, db):
    ln = open(file_name).read().replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '')

    read_block(ln, db)

    return db


class NodeTable:

    def __init__(self, label, sub_label="None"):
        self.label = label
        self.sub_label = sub_label
        self.data = []


def read_node_table(file_name, label, nodes=None):
    fin = open(file_name, 'r')

    start_label = str('<' + label)
    end_label = str('</' + label)

    output = []

    for line in fin:

        if line.strip().startswith(start_label):

            nt = NodeTable(label)

            if 'name' in line:
                sub_label = line.split('=')[1].replace('\n', '').replace('>', '').replace(' ', '').replace('\"', '').replace('\'', '')
                nt.sub_label = sub_label

            for line in fin:

                if line.strip().startswith(end_label):
                    output.append(nt)
                    break

                fullRel = line.strip().split(';')

                if len(fullRel) == 2:
                    splitRel = fullRel[0].split('=')

                    if len(splitRel) == 2:
                        lhs = splitRel[0]
                        rhs = splitRel[1]

                        if not is_node_dof(lhs):
                            raise RuntimeError(str(lhs) + ' is not a NodeDof')

                        dof_type, node_ids = decode_node_dof(lhs, nodes)

                        if get_type(rhs) is float or get_type(rhs) is int:
                            for node_id in node_ids:
                                nt.data.append([dof_type, int(node_id), float(eval(rhs))])
                        else:
                            rhs = rhs.replace(" ", "").replace("+", " +").replace("-", " -")
                            splitrhs = rhs.split(" ")
                            rhs = 0.0
                            for irhs in splitrhs:
                                if irhs == "":
                                    continue
                                if '[' not in irhs:
                                    for node_id in node_ids:
                                        nt.data.append([dof_type, int(node_id), float(eval(irhs))])
                                else:
                                    eq_rhs = irhs.split("*")
                                    factor = 1.0

                                    for ieq_rhs in eq_rhs:
                                        if (get_type(ieq_rhs) is float) or (get_type(ieq_rhs) is int):
                                            factor = clean_variable(ieq_rhs)
                                        else:
                                            if is_node_dof(ieq_rhs):
                                                if '-' in ieq_rhs:
                                                    factor = -1.0;
                                                ieq_rhs = ieq_rhs.replace("-", "").replace("+", "")
                                                slaveDofType, slaveNodeID = decode_node_dof(ieq_rhs, nodes)

                                    for node_id in node_ids:
                                        dt = [dof_type, int(node_id), rhs, slaveDofType, slaveNodeID, factor]
                                        nt.data.append(dt)

    return output
