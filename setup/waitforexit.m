function waitforexit( jobname, tag )
  doexit = 0;
  disp('Waiting...');
  xmlfilename = [tag, '.xml'];
  while ( doexit == 0 ) 
    system(['qstat -x ', jobname, ' > ', xmlfilename]);
    output = xml2struct( xmlfilename );
    for ii=1:length(output.Children.Children)
      if( strcmp( output.Children.Children(ii).Name, 'job_state') )
        if( strcmp(output.Children.Children(ii).Children.Data, 'C') )
          doexit = 1;
        end
      end
    end;
    pause(1);
  end;
  delete(xmlfilename);
  %  disp('Done.');
  
  