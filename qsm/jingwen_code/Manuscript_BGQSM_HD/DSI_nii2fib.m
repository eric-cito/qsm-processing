function [] = DSI_nii2fib(input_nii, output_fib)

addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));

nii_name = [input_nii '.nii.gz'];
fib_name = [output_fib '.fib'];

I = load_nii(nii_name);
image_size = size(I.img);
fib.dimension = image_size(1:3);
fib.voxel_size = I.hdr.dime.pixdim(2:4); % [1 1 1]; % [2.5 2.5 2.5]
fib.dir0 = zeros([3 fib.dimension]);
fib.fa0 = zeros(fib.dimension);
for z = 1:image_size(3)
    fprintf('# Slice %i\n', z);
    for y = 1:image_size(2)
        for x = 1:image_size(1)
            tensor = zeros(3,3);
            tensor(1,1) = I.img(x,y,z,1);
            tensor(1,2) = I.img(x,y,z,2);
            tensor(2,1) = I.img(x,y,z,2);
            tensor(1,3) = I.img(x,y,z,4);
            tensor(3,1) = I.img(x,y,z,4);
            tensor(2,2) = I.img(x,y,z,3);
            tensor(2,3) = I.img(x,y,z,5);
            tensor(3,2) = I.img(x,y,z,5);
            tensor(3,3) = I.img(x,y,z,6);
            [V D] = eig(tensor);
            if D(3,3) == 0
                continue;
            end
            l1 = D(3,3);
            if(l1 < 0)
                continue;
            end
            l2 = D(2,2);
            l3 = D(1,1);
            if(l2 < 0)
                l2 = 0;
            end
            if(l3 < 0)
                l3 = 0;
            end
            ll = (l1+l2+l3)/3;
            ll1 = l1-ll;
            ll2 = l2-ll;
            ll3 = l3-ll;
            fib.fa0(x,y,z) = sqrt(1.5*(ll1*ll1+ll2*ll2+ll3*ll3)/(l1*l1+l2*l2+l3*l3));
            V(1,3) = -V(1,3);
            fib.dir0(:,x,y,z) = V(:,3);
        end
    end
end

% enable the following codes if the image need to flip x and y
fib.fa0 = fib.fa0(image_size(1):-1:1,image_size(2):-1:1,:);
% fib.dir0 = fib.dir0(:,image_size(1):-1:1,image_size(2):-1:1,:);
% fib.dir0(3,:,:,:) = -fib.dir0(3,:,:,:);

fib.fa0 = reshape(fib.fa0,1,[]);
fib.dir0 = reshape(fib.dir0,3,[]);
save(fib_name,'-struct','fib','-v4');
gzip(fib_name);

end