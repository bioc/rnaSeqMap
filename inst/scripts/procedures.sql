DELIMITER $$

DROP PROCEDURE IF EXISTS nuc$$
CREATE PROCEDURE nuc(a int, b int)
BEGIN
declare k int default 0;
declare s int default 0;
set s = a;
set k = b;

DROP TEMPORARY TABLE IF EXISTS tmp ;  
CREATE TEMPORARY TABLE tmp ( n int(10));
insert into tmp values(s);

WHILE s < k DO
SET s=s+1;
insert into tmp values (s);
end while;

END$$


DROP PROCEDURE IF EXISTS seq_coverage$$
CREATE PROCEDURE seq_coverage(exp int, chr int, pstart int, pend int, pstr int)
BEGIN

call nuc(pstart,pend);

select tmp.n as "Nucleotide",count(sample_id) as "Number of reads" 
      from (select * from seq_read 
            where sample_id=exp
            and  chromosome_id=chr
            and  start >= pstart 
            and rend <= pend  
            and strand=pstr) as t right outer join tmp on( tmp.n between t.start and t.rend)  
      group by tmp.n;
END$$

DROP PROCEDURE IF EXISTS seq_histogram$$
CREATE PROCEDURE seq_histogram(exp int,chr int, s int, e int, str int, granular int)
BEGIN

declare a int default 0;
declare suma int default 0;
declare licz int default 0;
declare c int default 0;

DROP TEMPORARY TABLE IF EXISTS temp ;  
  CREATE TEMPORARY TABLE temp (
    n int(10), nr_reads int(10));
DROP TEMPORARY TABLE IF EXISTS t ;  
  CREATE TEMPORARY TABLE t (
    n int(10),c_reads int(10));

call nuc(s,e);

insert into temp
select tmp.n as "Nucleotide",count(*) as "Number of reads" 
      from (select * from seq_read 
            where sample_id=exp
            and chromosome_id=chr
            and  start >= s 
            and rend <= e 
            and strand=str) as t right outer join tmp on( tmp.n between t.start and t.rend)  
      group by tmp.n;

select count(*) from temp into suma;

while licz <= suma do
set licz=licz+granular;
SET @skip=c; SET @numrows=granular; 
PREPARE STMT FROM 'insert into t select tt.n, sum(tt.nr_reads) from ( SELECT n, nr_reads FROM temp LIMIT ?, ?) as tt';
EXECUTE STMT USING @skip, @numrows;
set c=licz-1;
end while;

select * from t;
END$$

DROP PROCEDURE IF EXISTS seq_indeks$$
CREATE PROCEDURE seq_indeks(exp1 int,exp2 int, chr int, s int, e int, str int)
BEGIN

declare const float(5,2) default 0;
declare sfc1 float(5,2) default 0;
declare sfc2 float(5,2) default 0;
declare len float(5,2) default 0;

DROP TEMPORARY TABLE IF EXISTS temp1 ;  
  CREATE TEMPORARY TABLE temp1 (sample_id int(2), n int(10), nr_reads int(5));
DROP TEMPORARY TABLE IF EXISTS temp2 ;  
  CREATE TEMPORARY TABLE temp2 (sample_id int(2), n int(10), nr_reads int(5));
DROP TEMPORARY TABLE IF EXISTS tabl ;  
  CREATE TEMPORARY TABLE tabl ( rule_id int,n int(10),c_reads int(10));
DROP TEMPORARY TABLE IF EXISTS indeks ;  
  CREATE TEMPORARY TABLE indeks (n int(10), ind int(10));

call nuc(s,e);

insert into temp1
select t.sample_id, tmp.n as "Nucleotide",log2(count(*)) as "nr_reads" 
      from (select * from seq_read 
            where sample_id = exp1 
            and chromosome_id = chr
            and  start >= s 
            and rend <= e 
            and strand=str) as t right outer join tmp on( tmp.n between t.start and t.rend)       
     group by tmp.n,t.sample_id;

insert into temp2
select t.sample_id, tmp.n as "Nucleotide",log2(count(*)) as "nr_reads" 
      from (select * from seq_read 
            where sample_id = exp2 
            and chromosome_id = chr
            and  start >= s 
            and rend <= e 
            and strand=str) as t right outer join tmp on( tmp.n between t.start and t.rend)       
     group by tmp.n,t.sample_id;

select count(*) from temp1 into len;
select sum(nr_reads) from temp1 into sfc1;
select sum(nr_reads) from temp2 into sfc2;

set const = (sfc1-sfc2)/len;
select n, truncate((temp1.nr_reads-temp2.nr_reads)/(const),3) from temp1 join temp2 using(n);
END$$

DROP PROCEDURE IF EXISTS readsInRange$$
CREATE PROCEDURE readsInRange(ex int, chr int, s int, e int, str int)
BEGIN

select start as "1",rend as "2" from seq_read 
where sample_id = ex
 and chromosome_id = chr
 and  start >= s 
 and rend <= e  
 and strand = str ;

END$$

DROP PROCEDURE IF EXISTS geneInChromosome$$
CREATE PROCEDURE geneInChromosome(chr varchar(5), s int, e int,str tinyint)
BEGIN

select symbol,stable_id,start,end,strand
from gene 
            where chromosome_name = chr
            and  start >= s 
            and end <= e 
            and strand = str
order by start;  

END$$


DROP PROCEDURE IF EXISTS t_gene$$
CREATE PROCEDURE t_gene(chr varchar(5), s int, e int,str tinyint)
BEGIN

DROP TEMPORARY TABLE IF EXISTS t_gene;  
CREATE TEMPORARY TABLE t_gene (stable_id varchar(128), symbol varchar(128), start int(10), stop int(10), strand tinyint(2));

insert into t_gene
select stable_id,symbol,start,end,strand
from gene 
            where chromosome_name = chr
            and  start >= s 
            and end <= e 
            and strand = str
order by start;  

END$$

DROP PROCEDURE IF EXISTS spaceInChromosome$$
CREATE PROCEDURE spaceInChromosome(chr varchar(5), s int, e int,str tinyint)
BEGIN

call t_gene(chr,s,e,str);
DROP TEMPORARY TABLE IF EXISTS t_space;  
  CREATE TEMPORARY TABLE t_space (start int(10), stop int(10), strand tinyint(2));

DROP TEMPORARY TABLE IF EXISTS t_space_p1;  
  CREATE TEMPORARY TABLE t_space_p1 (start int(10), stop int(10), znacznik int(10) PRIMARY KEY AUTO_INCREMENT);
DROP TEMPORARY TABLE IF EXISTS t_space_p2;  
  CREATE TEMPORARY TABLE t_space_p2 (start int(10), stop int(10), znacznik int(10) PRIMARY KEY AUTO_INCREMENT);
DROP TEMPORARY TABLE IF EXISTS tt;  
  CREATE TEMPORARY TABLE tt (start int(10), stop int(10), strand tinyint(2));
DROP TEMPORARY TABLE IF EXISTS ttt;  
  CREATE TEMPORARY TABLE ttt (start int(10), stop int(10), strand tinyint(2));

insert into tt select start,stop,strand from t_gene;
insert into ttt select start,stop,strand from t_gene;
insert into t_space_p2(start) values (s);

insert into t_space_p1(stop)
select distinct start-1 from t_gene where start-1 not in
(select tt.start-1 from tt, ttt where tt.start-1 between ttt.start and ttt.stop);

insert into t_space_p2(start)
select distinct stop+1 from t_gene where stop+1 not in
(select tt.stop+1 from tt, ttt where tt.stop+1 between ttt.start and ttt.stop);

insert into t_space_p1(stop) values (e);

insert into t_space
select p2.start,p1.stop,str from t_space_p1 p1 join t_space_p2 p2 using (znacznik);

select * from t_space;
drop temporary table tt;
drop temporary table ttt;
drop temporary table t_space_p1;
drop temporary table t_space_p2;
DROP TEMPORARY TABLE T_GENE;
END$$

DROP PROCEDURE IF EXISTS showDescription$$
CREATE  PROCEDURE showDescription()
BEGIN
select * from bio_sample;
END$$

DROP PROCEDURE IF EXISTS readsForGene$$
CREATE PROCEDURE readsForGene(chr int, st int, en int, str int)
BEGIN
 
DECLARE done INT DEFAULT 0 ;
DECLARE g_s char(20);
DECLARE g_start INT DEFAULT 0 ;
DECLARE g_end INT DEFAULT 0 ;
DECLARE g_strand INT DEFAULT 0 ;
 
DECLARE g_cur CURSOR FOR
  select stable_id, start, stop, strand from t_gene;

DECLARE CONTINUE HANDLER FOR NOT FOUND SET done = 1;
  
call t_gene(chr,st,en,str); 

DROP TEMPORARY TABLE IF EXISTS temp ;  
  CREATE TEMPORARY TABLE temp (symbol_id char(20), sample_id int(10), number_of_reads int(10));

OPEN g_cur; 

REPEAT      
FETCH g_cur INTO g_s, g_start, g_end, g_strand;           

IF NOT done THEN
 insert into temp
SELECT g_s, sample_id, count(sample_id) 
      from seq_read    
      where 
      chromosome_id=chr 
      and start >= g_start
      and rend <= g_end
      and strand = g_strand
group by sample_id;

END IF;

UNTIL done END REPEAT;

CLOSE g_cur;
    
select * from temp right outer join t_gene on (symbol_id=stable_id);

drop temporary table temp;
drop temporary table t_gene;
END$$


DROP PROCEDURE IF EXISTS seq_count_exps;  
CREATE PROCEDURE seq_count_exps()
BEGIN
select count(sample_id)
from seq_read
group by sample_id;

END$$

DELIMITER ;
