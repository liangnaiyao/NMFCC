function [U_final, V_final, objhistory_final] = obj_No_duibi(tryNo,U,V,objhistory,objhistory_final)
    if tryNo == 1       % ��һ�ε������뱻�Ƚ϶���
        U_final = U;
        V_final = V;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end) % �����ǵ�һ�ε��������ｲ�Ҽ���nRepeat��Ŀ�꺯����С���Ǵε����ķֽ���
           U_final = U;
           V_final = V;
           objhistory_final = objhistory;
       end
    end