import fractions

def mediant(a: fractions.Fraction, b: fractions.Fraction) -> fractions.Fraction:
    return fractions.Fraction(a.numerator + b.numerator, a.denominator + b.denominator)

def mediantSeries(n, m) -> fractions.Fraction:
    if n == 0:
        return fractions.Fraction(m, 1)
    else:
        w = m % 1
        # return (mediantSeries(n-1, m//2) + w * mediantSeries(n-1, m//2 + 1))
        if w:
            return mediant(mediantSeries(n-1, m//2), mediantSeries(n-1, m//2 + 1))
        else:
            return mediantSeries(n-1, m//2)

print(f'{mediantSeries(1, 1).numerator} / {mediantSeries(1,1).denominator}')


# breakpoint()
# for n in range(1):
#     for m in range(2 ** (n + 1) + 1):
#         x = mediantSeries(n,m)
#         print(f'{x.numerator} / {x.denominator}', end=', ')
#     print()