double precision function gaussian_sphere(x,y,z,w,a,hmax)
   implicit none

   double precision x,y,z,w(3),a,q, hmax

   double precision get_gc_distance, d

   d = get_gc_distance(x,y,z,w)
   q = hmax*exp(-a*d**2)

   gaussian_sphere = q

end

!! ------------------------------------------------
!! # Great circle distance on sphere
!! ------------------------------------------------
double precision function get_gc_distance(x,y,z,w)
   implicit none

   double precision x,y,z,w(3)
   double precision p(3), pn, wn, ca, th
   double precision rsphere, get_scale
   integer m

   rsphere = 1.d0

   p(1) = x
   p(2) = y
   p(3) = z
   pn = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
   if (abs(pn - rsphere)/rsphere .ge. 1e-3) then
      write(6,*) 'get_gc_distance : Point is not on the sphere; '
      write(6,*) x, y, z, abs(pn-rsphere)/rsphere
      stop
   endif

   wn = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
   do m = 1,3
      w(m) = w(m)/wn
   enddo

   ca = (w(1)*x + w(2)*y + w(3)*z)/pn
   if (abs(ca) .gt. 1) then
      write(6,*) 'get_gc_distance : abs(ca) > 1; ', ca
      write(6,*) w(1), w(2), w(3), x,y,z
      stop
   endif
   th = acos(ca)
   get_gc_distance = th*rsphere

end
