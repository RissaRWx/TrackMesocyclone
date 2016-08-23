load "./contributed.ncl_beta_640"

function track_meso(uh[*][*][*]:numeric,uh_thresh:numeric,minarea:numeric,dx:numeric,dt:numeric,storm_motion[*]:numeric)
local uh_high, ii_x, ii_y, nii_x, nii_y, dims, uhhi_1d_x, uhhi_1d_y, nmx, nmy, i, j, ibeg, iend, jbeg, jend, icn, jcn, unique_x, num_unique_x, unique_x_count, jj, dist_bw_x, m, num_objs, n, x_track, y_track, meso_tracks
begin
 
  dims = dimsizes(uh)
  times = dims(0)
  dim1 = dims(1)
  dim2 = dims(2)
  labels = new(dims,float)
  centers = new((/times,40,2/),float)
  centers_final = new((/times,40,2/),float)
  num_previous_tracks = 0  
  new_label = 1
do t = 0,times-1,1
  do y = 0,dim2-1,1
    do x = 0,dim1-1,1
      
      ; Go to the next point right if we don't meet the UH criteria 
      if(uh(t,y,x) .lt. uh_thresh) then
	continue
      else if (y .gt. 0 .and. uh(t,y-1,x) .ge. uh_thresh) then
	labels(t,y,x) = labels(t,y-1,x)
      else if (x+1 .lt. dim2 .and. y .gt. 0 .and. uh(t,y-1,x+1) .ge. uh_thresh) then
	c = labels(t,y-1,x+1)
	labels(t,y,x) = c
	if (x .gt. 0 .and. uh(t,y-1,x-1) .ge. uh_thresh) then
	  a = labels(t,y-1,x-1)
	  if (c .ne. a) then
	    labels(t,:,:) = where(labels(t,:,:) .eq. c, a, labels(t,:,:))
	  end if 
	else if (x .gt. 0 .and. uh(t,y,x-1) .ge. uh_thresh) then
	  d = labels(t,y,x-1)
	  if (d .ne. c) then
	    labels(t,:,:) = where(labels(t,:,:) .eq. d, c, labels(t,:,:))
	  end if
        end if
	end if
      else if (x .gt. 0 .and. y .gt. 0 .and. uh(t,y-1,x-1) .ge. uh_thresh) then
	labels(t,y,x) = labels(t,y-1,x-1)
      else if (x .gt. 0 .and. uh(t,y,x-1) .ge. uh_thresh) then
	labels(t,y,x) = labels(t,y,x-1)
      else
	labels(t,y,x) = new_label
	new_label = new_label + 1
      end if
      end if
      end if
      end if
      end if
    end do
  end do
  labels_list = get_unique_values(labels(t,:,:))
  if (dimsizes(labels_list) .eq. 1) then
    delete(labels_list)
    continue
  else
    num_labels = count_unique_values(labels(t,:,:))
    label_counts = new(num_labels,float)
    do i = 0,num_labels-1
      label_counts(i) = num(labels(t,:,:).eq.labels_list(i))
    end do
    small_areas_labels = labels_list(ind(label_counts .lt. minarea))
    delete([/labels_list,label_counts/])
    if (.not.all(ismissing(small_areas_labels)))
      do i = 0,dimsizes(small_areas_labels)-1,1
        labels(t,:,:) = where(labels(t,:,:) .eq. small_areas_labels(i), labels@_FillValue,labels(t,:,:))
      end do
      delete(small_areas_labels)
    end if
    labels_list_new = get_unique_values(labels(t,:,:))
    if (all(ismissing(labels(t,:,:)))) then		; If we haven't got any areas larger than min area, go to the next time
      delete(labels_list_new)
      continue
    end if
    do i = 0,dimsizes(labels_list_new)-1,1
      labels_1d = ndtooned(labels(t,:,:))
      locs_1d = ind(labels_1d .eq. labels_list_new(i))
      locs = ind_resolve(locs_1d,dimsizes(labels(t,:,:)))
      centers(t,i,0) = ceil((max(locs(:,0))+min(locs(:,0)))/2)
      centers(t,i,1) = ceil((max(locs(:,1))+min(locs(:,1)))/2)
      delete([/locs,locs_1d/])
    end do
    delete(labels_list_new)
    if(t.eq. 0 .or. all(ismissing(centers_final(t-1,:,:))))
      num_new_centers = num(.not.ismissing(centers(t,:,0)))
      centers_final(t,num_previous_tracks:toint(num_previous_tracks+num_new_centers-1),:) = centers(t,0:toint(num_new_centers-1),:)
      num_previous_tracks = num_previous_tracks + num_new_centers

      continue
    else
      xdist = ceil(storm_motion(0)*dt/dx)
      ydist = ceil(storm_motion(1)*dt/dx)*-1
      centers_pred_y = centers_final(t-1,:,0)+ydist
      centers_pred_x = centers_final(t-1,:,1)+xdist
      num_centers_current = num(.not.ismissing(centers(t,:,0)))
      current_centers_ind = ind(.not.ismissing(centers(t,:,0)))
      do i = 0,num(.not.ismissing(centers(t,:,0)))-1
         good_center_cur = current_centers_ind(i)
        ;if (num_centers_current .gt. 1) then
          diff = sqrt((centers(t,good_center_cur,0)*dx-centers_pred_y(:)*dx)^2+(centers(t,good_center_cur,1)*dx-centers_pred_x(:)*dx)^2)
          min_diff = min(diff)
          min_diff_center_ind = ind(diff .eq. min_diff)
          min_diff_center = min_diff_center_ind(0)
          delete(min_diff_center_ind)
        if (min_diff .lt. 7500.0) then	; We're within 5km of the predicted track, consider associating with a previous track location
          if (.not.ismissing(centers_final(t,min_diff_center,0))) then	; If we've already got a center on the nearest track at this time, see whether old or new center is closer
            min_diff_old = sqrt((centers_final(t,min_diff_center,0)*dx-centers_pred_y(min_diff_center)*dx)^2+(centers_final(t,min_diff_center,1)*dx-centers_pred_x(min_diff_center)*dx)^2)
            min_diff = min((/min_diff,min_diff_old/))
            if (min_diff .eq.  min_diff_old) then	; If old center is closer, ignore the new center to prevent tracks from being < 5km from each other
              continue
            end if 
          end if
          centers_final(t,min_diff_center,:) = centers(t,good_center_cur,:)
        else
          centers_final(t,num_previous_tracks,:) = centers(t,good_center_cur,:)
          num_previous_tracks = num_previous_tracks + 1
        end if
      end do
      delete([/current_centers_ind,min_diff_center,centers_pred_y,centers_pred_x/])
    end if
  end if
end do
return(centers_final)
end  


