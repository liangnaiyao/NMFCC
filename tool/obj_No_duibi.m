function [U_final, V_final, objhistory_final] = obj_No_duibi(tryNo,U,V,objhistory,objhistory_final)
    if tryNo == 1       % 第一次迭代存入被比较对象
        U_final = U;
        V_final = V;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end) % 当不是第一次迭代，这里讲找几次nRepeat中目标函数最小的那次迭代的分解结果
           U_final = U;
           V_final = V;
           objhistory_final = objhistory;
       end
    end