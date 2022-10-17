//namespace WignerSymbols {

template <typename T>
double sgn(T val)
{
    int sgn = (T(0) < val) - (val < T(0));
    if (sgn == 0)
        return 1.0;
    else
        return (double)sgn;
}

//} // namespace WignerSymbols
