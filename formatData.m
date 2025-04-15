function [matrix_pos_1, matrix_pos_2, matrix_neg_1, matrix_neg_2] = formatData(data)
    matrix_pos_1 = horzcat(data, zeros(size(data,1),1));
    matrix_pos_2 = horzcat(zeros(size(data,1),1), data);
    matrix_neg_1 = horzcat(-data, zeros(size(data,1),1));
    matrix_neg_2 = horzcat(zeros(size(data,1),1), -data);
end


