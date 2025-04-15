function [matrix_pos_1, matrix_pos_2,matrix_pos_3,matrix_neg_1,matrix_neg_2,matrix_neg_3] = formatData(data)
    matrix_pos_1 = horzcat(data, zeros(size(data,1),2));
    matrix_pos_2 = horzcat(zeros(size(data,1),1),data,zeros(size(data,1),1));
    matrix_pos_3 = horzcat(zeros(size(data,1),2), data);
    matrix_neg_1 = horzcat(-data, zeros(size(data,1),2));
    matrix_neg_2 = horzcat(zeros(size(data,1),1),-data, zeros(size(data,1),1));
    matrix_neg_3 = horzcat(zeros(size(data,1),2),-data);
end


