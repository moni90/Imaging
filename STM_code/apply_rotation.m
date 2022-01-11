function f_rotated = apply_rotation(f_pattern, teta)
%build rotation matrix
if teta == pi/2
    f_rotated = f_pattern;
else
    e0 = [1 0]';
    e1 = [cos(teta) sin(teta)]';
    A = [e0 e1];
    cos_t = cos(teta);
    sin_t = sin(teta);

    while size(A,1)<size(f_pattern,1)
        A = [A; zeros(1,size(A,1))];
        cos_t = (cos_t-cos_t^2)/(sin_t^2);
        if isnan(cos_t)
            cos_t = 0;
        end
        sin_t = sqrt(1-cos_t^2);
        e_new = A(:,end);
        e_new(end-1) = A(end-1,end)*cos_t;
        e_new(end) = A(end-1,end)*sin_t;
        A = [A e_new];
    end
    %apply rotation matrix to signals
    f_rotated = A*f_pattern;

end
