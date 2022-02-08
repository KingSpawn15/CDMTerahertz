
sender = 'rustomji@post.bgu.ac.il';
psswd = 'sg3juPms';


setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
filename = mfilename('path');
% sendmail('kaizad.rustomji@gmail.com','Finished','Hello! This is a test from MATLAB!')