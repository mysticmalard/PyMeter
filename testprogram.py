breakpoint()
a = 1

def foo(x, y):
    return (x, y)

print(globals())
# print(__annotations__)
print(__builtins__)
# print(__dict__)
print(type(__builtins__))
print(a.__dict__)