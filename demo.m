clear all
warning('off')


addpath('./qtfm/') % 

O_Img = double(imread('./starfish.png'));

nSig  = 40;

randn('seed', 0);
N_Img = O_Img + nSig* randn(size(O_Img));                                   %Generate noisy image
PSNR  =  csnr( N_Img, O_Img, 0, 0 );
fprintf( 'Noisy Image: nSig = %2.3f, PSNR = %2.2f \n\n\n', nSig, PSNR );


figure
Par   = Q_ParSet(nSig);   
E_Img = QNMF_DeNoising( N_Img, O_Img, Par );                                %WNNM denoisng function
E_Img = uint8(E_Img);

imshow(E_Img)
imwrite(E_Img , 'denoise.png')