% Produce bode plots of TCV plant
load TCV_model

axm = [];
axp = [];

out_names = output_names;
out_names(find(out_names == '_')) = ' ';

w = 2*pi*logspace(-3, 4);
G = freqresp(ss(A,B,C,D),w);
M = abs(G);
P = angle(G);

OH_set = [1:2];
E_set_a = [3:6];
E_set_b = [7:10];
F_set_a = [11:14];
F_set_b = [15:18];
G_set = [19];

for n=1:5
   if n==1 
      coil_set = OH_set;
   elseif n==2
      coil_set = E_set_a;
   elseif n==3
      coil_set = F_set_a;
   elseif n==4
      coil_set = E_set_b;
   elseif n==5
      coil_set = F_set_b;
   end
   figure(n)
   clf
   set(gca,'fontsize',8)
   for ins = coil_set;
      for outs = 1:6;
         %panelopen(length(coil_set)*2,find(ins==coil_set)*2-1, size(C,1),outs);
         subplot(length(coil_set),6, (outs + (find(ins==coil_set)-1)*6))
         m = M(outs,ins,:);m=m(:);
         semilogx(w, 20*log10(m), 'b')
         a = axis;
         axm = [axm;a];
         axis([min(w) max(w) -200 100])
         set(gca,'xtick',[1e-2 1e0 1e2 1e4])
         if ins == min(coil_set)
            title(out_names(outs,:))
         end
         if outs==1
            ylabel(coil_names(ins,:))            
         end
         if ins == max(coil_set)
            xlabel('\omega / rad s^{-1}')
         end
         grid on
         %panelclose(length(coil_set)*2,find(ins==coil_set)*2-1, size(C,1),outs);
         
         %panelopen(length(coil_set)*2,find(ins==coil_set)*2, size(C,1),outs);
         p = P(outs,ins,:);p=p(:);
         %semilogx(w, p, 'r')
         a = axis;
         axp = [axp;a];
         %axis([min(w) max(w) a(3) a(4)])
         %panelclose(length(coil_set)*2,find(ins==coil_set)*2, size(C,1),outs);
         
      end
   end
   orient landscape 
   %eval(['print -depsc2 D:/matlab/MPC/TCV/TCV_freqresp_' num2str(n) '.eps'])
end


