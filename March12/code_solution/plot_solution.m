cd output/
xaxis
bed

prompt = {'Give the number of the files'};
titolo = '';
dati = inputdlg(prompt,titolo,1);
name = char(dati(1));
nfiles=str2double(char(dati(1)));

for i=1:nfiles
    if(i<10)
        eval(['depth0' num2str(i)]);
    else
        eval(['depth0' num2str(i)]);
    end
    plot(x,b,x,h);
    pause(0.1)
end
depth
plot(x,b,x,h)
