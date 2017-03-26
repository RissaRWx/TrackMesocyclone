load "./contributed.ncl_beta_640"

function track_meso(uh[*][*][*]:numeric,vars,uh_thresh:numeric,minarea:numeric,dx:numeric,dt:numeric,storm_motion[*]:numeric,track_tol:numeric)
local uh_high, ii_x, ii_y, nii_x, nii_y, dims, uhhi_1d_x, uhhi_1d_y, nmx, nmy, i, j, ibeg, iend, jbeg, jend, icn, jcn, unique_x, num_unique_x, unique_x_count, jj, dist_bw_x, m, num_objs, n, x_track, y_track, meso_tracks,avo_meso_cur,avo_cur
begin
 
  dims = dimsizes(uh)
  times = dims(0)
  dim1 = dims(1)
  dim2 = dims(2)
  printVarSummary(vars)
  numvars = ListCount(vars)
  labels = new(dims,float)
  centers = new((/times,60,2/),float)
  centers_final = new((/times,60,2/),float)
  vars_meso = new((/numvars,times,60/),float)
  accum_vars = new((/numvars,times,60/),float)
  num_previous_tracks = 0  
  new_label = 1
  
do t = 0,times-1,1
  do y = 0,dim1-1,1
    do x = 0,dim2-1,1
      
      if(uh(t,y,x) .lt. uh_thresh) then  ; Go to the next point right if we don't meet the UH criteria
	continue
      else if (y .gt. 0 .and. uh(t,y-1,x) .ge. uh_thresh) then ; If B is gt threshold, give this point same label
	labels(t,y,x) = labels(t,y-1,x)
      else if (x+1 .lt. dim2 .and. y .gt. 0 .and. uh(t,y-1,x+1) .ge. uh_thresh) then ; If BR is labeled, must check other areas for labels
	c = labels(t,y-1,x+1)
	labels(t,y,x) = c
	if (x .gt. 0 .and. uh(t,y-1,x-1) .ge. uh_thresh) then ; If point to BL is labeled as well, BL and BR areas connect through this point
	  a = labels(t,y-1,x-1)
	  if (c .ne. a) then
	    labels(t,:,:) = where(labels(t,:,:) .eq. c, a, labels(t,:,:))  ; If BL and BR points are different areas, are now the same area
	  end if 
	else if (x .gt. 0 .and. uh(t,y,x-1) .ge. uh_thresh) then ; If L is labeled but not BL, check to see if L=BR
	  d = labels(t,y,x-1)
	  if (d .ne. c) then
	    labels(t,:,:) = where(labels(t,:,:) .eq. d, c, labels(t,:,:)) ; If B = BR, B and BR areas connect through P
	  end if
        end if
	end if
      else if (x .gt. 0 .and. y .gt. 0 .and. uh(t,y-1,x-1) .ge. uh_thresh) then ; If BL is labeled but BR and B aren't, P=BL
	labels(t,y,x) = labels(t,y-1,x-1)
      else if (x .gt. 0 .and. uh(t,y,x-1) .ge. uh_thresh) then ; If L is labeled but BR, B, and BL aren't, P=L
	labels(t,y,x) = labels(t,y,x-1)
      else		; BR, BL, B, and L aren't labeled, we've got a new region! Create a new label
	labels(t,y,x) = new_label
	new_label = new_label + 1
      end if
      end if
      end if
      end if
      end if
    end do
  end do
  labels_list = get_unique_values(labels(t,:,:)) ; Find how many unique region labels we have
  if (dimsizes(labels_list) .eq. 1) then ; If we don't have any regions, go to the next time step
    delete(labels_list)
    continue
  else
    num_labels = count_unique_values(labels(t,:,:)) ; Count how many unique labels we have
    label_counts = new(num_labels,float) ; 
    do i = 0,num_labels-1 ; Count how many times each label appears to determine the size of the area
      label_counts(i) = num(labels(t,:,:).eq.labels_list(i))
    end do
    if (.not.all(ismissing(ind(label_counts .lt. minarea))))
      small_areas_labels = labels_list(ind(label_counts .lt. minarea)) ; Find areas smaller than the threshold size
      delete([/labels_list,label_counts/])
      if (.not.all(ismissing(small_areas_labels)))
        do i = 0,dimsizes(small_areas_labels)-1,1
          labels(t,:,:) = where(labels(t,:,:) .eq. small_areas_labels(i), labels@_FillValue,labels(t,:,:)) ; Remove labels for small areas
        end do
        delete(small_areas_labels)
      end if
    else delete([/labels_list,label_counts/])
    end if
    labels_list_new = get_unique_values(labels(t,:,:))
    if (all(ismissing(labels(t,:,:)))) then		; If we haven't got any areas larger than min area, go to the next time
      delete(labels_list_new)
      continue
    end if
    uh_current = ndtooned(uh(t,:,:))
    do i = 0,dimsizes(labels_list_new)-1,1
      labels_1d = ndtooned(labels(t,:,:))
      locs_1d = ind(labels_1d .eq. labels_list_new(i))
      locs = ind_resolve(locs_1d,dimsizes(labels(t,:,:))) ; find the x-y pairs for the indices with the current label
      centers(t,i,0) = ceil(sum(locs(:,0)*uh_current(locs_1d))/sum(uh_current(locs_1d)))
      centers(t,i,1) = ceil(sum(locs(:,1)*uh_current(locs_1d))/sum(uh_current(locs_1d)))
      do j=0,numvars-1
        var_cur = vars[j]
        var_cur := ndtooned(var_cur(t,:,:))
        vars_meso(j,t,i) = avg(var_cur(locs_1d))
        delete(var_cur)
      end do
      ;centers(t,i,0) = ceil((max(locs(:,0))+min(locs(:,0)))/2)
      ;centers(t,i,1) = ceil((max(locs(:,1))+min(locs(:,1)))/2)
      delete([/locs,locs_1d/])
    end do
    delete(labels_list_new)
    if(t.eq. 0 .or. all(ismissing(centers_final(t-1,:,:)))) ; if we're at the first time step or there are no centers at the previous time, all IDd centers are new tracks
      num_new_centers = num(.not.ismissing(centers(t,:,0)))
      centers_final(t,num_previous_tracks:toint(num_previous_tracks+num_new_centers-1),:) = centers(t,0:toint(num_new_centers-1),:)
      num_previous_tracks = num_previous_tracks + num_new_centers

      continue
    else
      xdist = ceil(storm_motion(0)*dt/dx)	; Dist traveled in x-dir based on storm motion
      ydist = ceil(storm_motion(1)*dt/dx)*-1	; Dist traveled in y-dir based on storm motion
      centers_pred_y = centers_final(t-1,:,0)+ydist	; Guess new center locs w/storm motion mod
      centers_pred_x = centers_final(t-1,:,1)+xdist
      num_centers_current = num(.not.ismissing(centers(t,:,0)))		; # of centers found at current time
      current_centers_ind = ind(.not.ismissing(centers(t,:,0)))		; Location of current centers in array
      do i = 0,num(.not.ismissing(centers(t,:,0)))-1
         good_center_cur = current_centers_ind(i)
        ;if (num_centers_current .gt. 1) then
          diff = sqrt((centers(t,good_center_cur,0)*dx-centers_pred_y(:)*dx)^2+(centers(t,good_center_cur,1)*dx-centers_pred_x(:)*dx)^2)	; Distance from current to predicted center
          min_diff = min(diff)
          min_diff_center_ind = ind(diff .eq. min_diff)
          min_diff_center = min_diff_center_ind(0)
          delete(min_diff_center_ind)
        ;else 
        ;  diff = sqrt((centers(t,i,0)*dx-centers_pred_y*dx)^2+(centers(t,i,1)*dx-centers_pred_x*dx)^2)
        ;  min_diff_center = ind(.not.ismissing(centers_final(t-1,:,0)))
        ;end if
        if (min_diff .lt. track_tol) then	; We're within track_tol of the predicted track, consider associating with a previous track location
          if (.not.ismissing(centers_final(t,min_diff_center,0))) then	; If we've already got a center on the nearest track at this time, see whether old or new center is closer
            min_diff_old = sqrt((centers_final(t,min_diff_center,0)*dx-centers_pred_y(min_diff_center)*dx)^2+(centers_final(t,min_diff_center,1)*dx-centers_pred_x(min_diff_center)*dx)^2)
            min_diff = min((/min_diff,min_diff_old/))
            if (min_diff .eq.  min_diff_old) then	; If old center is closer, ignore the new center to prevent tracks from being < 5km from each other
              continue
            end if 
          end if
          centers_final(t,min_diff_center,:) = centers(t,good_center_cur,:) ; Set final center location
          accum_vars(:,t,min_diff_center) = vars_meso(:,t,good_center_cur)* \
					sqrt((centers_final(t,min_diff_center,0)*dx-centers_final(t-1,min_diff_center,0)*dx)^2+ \
					(centers_final(t,min_diff_center,1)*dx-centers_final(t-1,min_diff_center,1)*dx)^2)
        else
          centers_final(t,num_previous_tracks,:) = centers(t,good_center_cur,:) ; Create a new track and set final center location
          accum_vars(:,t,num_previous_tracks) = 0.0
          num_previous_tracks = num_previous_tracks + 1
        end if
      end do
      delete([/current_centers_ind,min_diff_center,centers_pred_y,centers_pred_x/])
    end if
  end if
end do
return([/centers_final,vars_meso,accum_vars/])
end  


