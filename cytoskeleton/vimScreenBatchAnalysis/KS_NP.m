 function ks_np = KS_NP(ks)
        
        ks_flip = ks';
      ks_np = max(ks, ks_flip);
      ks_np(find(ks<ks_flip)) = - ks_np(find(ks<ks_flip));
