import operator
import math
import copy
import matplotlib.pyplot as plt

OPERATION_PREDECENCE = {"*": 3, "^": 4, "/": 3, "+": 2, "-": 2, "sin": 1, "max": 1, "min": 1, "cos": 1, "tan": 1, "sqrt": 1, "abs": 1, "floor": 1, "fact": 1}
OPERATOR_DIC = {"+": operator.add, "^": operator.pow, "/": operator.truediv, "*": operator.mul, "-": operator.sub}
ACCEPTED_CARACTERS = "+/*-()^"
FUNCTIONS = {"sin": math.sin, "cos": math.cos, "tan": math.tan, "sqrt": math.sqrt, "abs": abs, "max": max, "min": min, "floor": math.floor, "fact": math.factorial}
VARNUMBER =  {"sin": 1, "cos":1, "tan": 1, "sqrt": 1, "abs": 1, "max": 2, "min": 2, "floor": 1, "fact": 1}
CONSTANTS = {"e": math.e, "pi": math.pi}

def arrange(lowerBound, upperBound, steps):
    lst = []
    cnt = lowerBound
    while cnt <= upperBound:
        lst.append(cnt)
        cnt += steps

    return lst
    

def isfloat(element):
    try :
        float(element)
        return True
    except:
        return False

def isint(element):
    try :
        int(element)
        return True
    except:
        return False
    
def parser_expression(expression, symbols):
    """
    :param expression: str,
    :param symbols: list,
    :return: The parsed expression as a list of operation, symbols and numbers for the Shuting Yard Algorithm
    """
    
    execution_task = []
    current_task = ""
    for i in expression:
        if i != " ":
            if i.isalnum() or i ==".":
                current_task += i
            elif i == ",":
                if isint(current_task):
                    execution_task.append(int(current_task))
                elif isfloat(current_task):
                    execution_task.append(float(current_task))
                current_task = ""
            elif i in ACCEPTED_CARACTERS:
                if current_task != "":
                    if isint(current_task):
                        execution_task.append(int(current_task))
                    elif isfloat(current_task):
                        execution_task.append(float(current_task))
                    elif (current_task in symbols) or (current_task in FUNCTIONS) or (current_task in CONSTANTS):
                        execution_task.append(current_task)
                    else:
                        raise ValueError(f"One of the symbols doesn't appear in the symbols list : {current_task}")
                execution_task.append(i)
                current_task = ""
    if current_task != "":
        if isint(current_task):
            execution_task.append(int(current_task))
        elif isfloat(current_task):
            execution_task.append(float(current_task))
        else:
            execution_task.append(current_task)
    return execution_task

def shuting_yard(expression):
    stack = []
    operator_stack = []
    while expression != []:
        i = expression.pop(0)
        if (type(i) in (float, int)) or (i not in ACCEPTED_CARACTERS and i not in FUNCTIONS) or (i in CONSTANTS):
            stack.append(i)
        elif i in FUNCTIONS:
            operator_stack.append(i)
        elif i in "*+^/-":
            if operator_stack == []: operator_stack.append(i)
            else:
                while operator_stack != [] and operator_stack[-1] != "(" and OPERATION_PREDECENCE[i] <= OPERATION_PREDECENCE[operator_stack[-1]]:
                    if i == "^" and operator_stack[-1] == "^":
                        break
                    else:
                        stack.append(operator_stack.pop())
                operator_stack.append(i)
        elif i == "(":
            operator_stack.append(i)
        elif i == ")":
            if operator_stack == []: raise SyntaxError("Mismatching parenthesis : '(' parenthesis missing")
            while operator_stack[-1] != "(":
                stack.append(operator_stack.pop())
                if operator_stack == [] : raise SyntaxError("Mismatching parenthesis : '(' parenthesis missing")
            operator_stack.pop()
            if (operator_stack != []) and (operator_stack[-1] in FUNCTIONS):
                stack.append(operator_stack.pop())

    while operator_stack != []:
        if operator_stack[-1] == "(": raise SyntaxError("Mismatching parenthesis : ')' parenthesis missing")
        stack.append(operator_stack.pop())
    return stack


def evaluate(rpn, symbols):
    exp = copy.deepcopy(rpn)
    stack = []
    while exp != []:
        i = exp.pop(0)
        if type(i) in (float, int):
            stack.append(i)
        elif i in FUNCTIONS:
            if VARNUMBER[i] == 1:
                stack.append(FUNCTIONS[i](stack.pop()))
            else:
                var = FUNCTIONS[i]([stack.pop() for _ in range(VARNUMBER[i])])
                stack.append(var)
        elif i in ACCEPTED_CARACTERS:
                right = stack.pop()
                left = stack.pop()
                stack.append(OPERATOR_DIC[i](left, right))
        elif i in CONSTANTS:
            stack.append(CONSTANTS[i])
        elif i in symbols:
            stack.append(symbols[i])

    return stack.pop()


class MathFunctions():
    def __init__(self, expression, symbols):
        self.symbols = symbols
        self.StringExpression = expression
        self.ParsedExpression = parser_expression(expression, symbols)
        self.ReadyToEvaluateExp = shuting_yard(self.ParsedExpression)

    def __str__(self):
        return self.StringExpression + ", Symbols : " + str(self.symbols)

    def eval(self, symbols_values):
        return evaluate(self.ReadyToEvaluateExp, symbols_values)

    def drawInRange2D(self, lowerBound, upperBound, steps):
        fig, axs = plt.subplots()
        ax = []
        ay = []
        for i in arrange(lowerBound, upperBound, steps):
            ax.append(i)
            ay.append(self.eval({self.symbols[0]: i}))
        axs.plot(ax, ay)
        plt.title(str(self))
        plt.grid(True)
        plt.show()
"""
x2 = MathFunctions("3*x^3 + 2*x - 1", ["x"])

x2.drawInRange2D(-5, 5, 0.001)"""

exp = MathFunctions(input("function: "), input("symbols: ").split())
exp.drawInRange2D(float(input("lower bound: ")), float(input("upper bound: ")), float(input("steps: ")))
