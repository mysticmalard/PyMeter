# import oldpymeter as pymeter
import pymeter1 as pymeter

m = pymeter.Value(42, 'kg')

c = pymeter.Value(299792458, 'm / s')

# breakpoint()
E = m * c ** 2
# breakpoint()
print(E)

# breakpoint()
d = pymeter.Value(10, 'm')

t = pymeter.Value(2, 's')

v = d / t

print(v)

d = pymeter.Value(13.6, 'm')

F = pymeter.Value(13000, 'N')

print(F * d)