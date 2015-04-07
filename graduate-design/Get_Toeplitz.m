function G = Get_Toeplitz(B, img_size, ker_size)

G = zeros((img_size(1)-ker_size(1)+1)*(img_size(2)-ker_size(2)+1), ker_size(1)* ker_size(2)); 
w = ker_size(1);
h = ker_size(2);

for i = 1:img_size(2)-ker_size(2)+1
    for j = 1:img_size(1)-ker_size(1)+1
        sub_mat = B(j+w-1:-1:j, i+h-1:-1:i);
        G((img_size(1)-ker_size(1) +1)*(i-1) + j, :) = sub_mat(:)';
    end
end

end

