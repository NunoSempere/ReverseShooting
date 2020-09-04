% nberdatesQ.m  1/28/09  -- by Quarter rather than Month
%  Official NBER recession datesQ, from
%   http://www.nber.org/cycles/cyclesmain.html  (see NBERcycles.html)

datesQ=[
% Peak Year,Quarter   Trough Year, Quarter
1857 2     1858 4;
1860 3     1861 3;
1865 1     1867 1;
1869 2     1870 4;
1873 3     1879 1;
1882 1     1885 2;
1887 2     1888 1;
1890 3     1891 2;
1893 1     1894 2;
1895 4     1897 2;
1899 3     1900 4;
1902 4     1904 3;
1907 2     1908 2;
1910 1     1912 4;
1913 1     1914 4;
1918 3     1919 1;
1920 1     1921 3;
1923 2     1924 3;
1926 3     1927 4;
1929 3     1933 1;
1937 2     1938 2;
1945 1     1945 4;
1948 4     1949 4;
1953 2     1954 2;
1957 3     1958 2;
1960 2     1961 1;
1969 4     1970 4;
1973 4     1975 1;
1980 1     1980 3;
1981 3     1982 4;
1990 3     1991 1;
2001 1     2001 4;
2007 4     2009 2       ];

start=datesQ(:,1)+(datesQ(:,2)-1)/4;
finish=datesQ(:,3)+(datesQ(:,4)-1)/4;
save nberdatesQ;































