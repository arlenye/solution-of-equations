from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30

#参数化解集
class Parametrization(object):
    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM_MSG = (
        'The basepoint and direction vectors should all live in the same dimension')

    def __init__(self,basepoint,direction_vectors):
        self.basepoint=basepoint
        self.direction_vectors=direction_vectors
        self.dimension=self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension==self.dimension

        except AssertionError:
            raise Exception(BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM_MSG)

class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    #交换两行
    def swap_rows(self, row1, row2):
        self.planes[row1],self.planes[row2] = self.planes[row2],self.planes[row1]

    #将等式乘以非零数字
    def multiply_coefficient_and_row(self, coefficient, row):
        n = self[row].normal_vector
        k = self[row].constant_term

        new_normal_vector = n.times_scalar(coefficient)
        new_constant_term = k * coefficient

        self[row] = Plane(normal_vector = new_normal_vector,constant_term = new_constant_term)

    #将多倍的等式加到另一个等式上
    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = (k1 * coefficient) + k2

        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector,constant_term=new_constant_term)

    #负责找出每个等式的第一个非零项
    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    #返回的是方程组里的平面数量
    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    #输出整洁版本的方程组
    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    #三角化方程组
    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension

        j=0
        for i in range(num_equations):
            while j<num_variables:
                c=MyDecimal(system[i].normal_vector.coordinates[j])
                if c.is_near_zero():
                    swap_succeeded=system.swap_with_row_below_for_nonzero_coefficient_if_able(i,j)
                    if not swap_succeeded:
                        j+=1
                        continue
                system.clear_coefficients_below(i,j)
                j+=1
                break

        return system

    #若第row行下面有某一行的j变量不为零，则交换这两行
    def swap_with_row_below_for_nonzero_coefficient_if_able(self,row,col):
        num_equations = len(self)

        for k in range(row+1,num_equations):
            coefficient = MyDecimal(self[k].normal_vector[col])
            if not coefficient.is_near_zaro():
                self.swap_rows(row,k)
                return True

        return False


    def clear_coefficients_below(self,row,col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector.coordinates[col])

        for k in range(row+1,num_equations):
            n = self[k].normal_vector
            gamma = n.coordinates[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha,row,k)

    #将方程组变为最简化的梯阵形式
    def compute_rref(self):
        tf = self.compute_triangular_form()

        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(num_equations)[::-1]:
            j=pivot_indices[i]
            if j<0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(i,j)
            tf.clear_coefficients_above(i,j)

        return tf

    #化系数为1
    def scale_row_to_make_coefficient_equal_one(self,row,col):
        n = self[row].normal_vector
        beta = Decimal('1.0')/n.coordinates[col]
        self.multiply_coefficient_and_row(beta,row)

    #从下往上消除系数
    def clear_coefficients_above(self,row,col):
        for k in range(row)[::1]:
            n=self[k].normal_vector
            alpha = -(n.coordinates[col])
            self.add_multiple_times_row_to_row(alpha,row,k)

    #求方程组的解
    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_parametrize_solution()

        except Exception as e:
            if str(e)==self.NO_SOLUTIONS_MSG:
                return str(e)
            else:
                raise e
            
    #找出唯一解
    def do_gaussian_elimination_and_parametrize_solution(self):
        rref = self.compute_rref()

        rref.raise_exception_if_contradictory_equation()

        direction_vectors = rref.extract_direction_vectors_for_parametrization()
        basepoint = rref.extract_basepoint_for_parametrization()

        return Parametrization(basepoint,direction_vectors)

    #求得参数化形式的方向向量
    def extract_direction_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        #free_variable_indices = set(range(num_variables)) - set(pivot_indices)
        
        i = 0
        all_variables = []
        while i<num_variables:
            all_variables.append(i)
            i+=1
        free_variable_indices = []
        for e in all_variables:
            if e not in pivot_indices:
                free_variable_indices.append(e)
                
        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0]*num_variables
            vector_coords[free_var] = 1
            for i,p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if pivot_var<0:
                    break
                vector_coords[pivot_var] = -p.normal_vector.coordinates[free_var]
            direction_vectors.append(Vector(vector_coords))
        
        return direction_vectors

    #求得基准点
    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0]*num_variables

        for i,p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if pivot_var<0:
                break
            basepoint_coords[pivot_var] = p.constant_term

        return Vector(basepoint_coords)

    #检查是否存在有矛盾的等式。检查每个平面，看其法向量的坐标是否全为0，如果全为0，则寻找常量项的非零项
    #如果能够找到，则说明存在矛盾的等式(即0=k)
    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)

            except Exception as e:
                if str(e)=='No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    #检查是否具有太多主变量
    def raise_exception_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index>=0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots<num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

        
            

        

#快速检查小数对象是否位于零公差范围内
class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

